package task

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/Taichidasheen/site_predict/pkg/calc"
	"github.com/Taichidasheen/site_predict/pkg/predict"
	tf "github.com/wamuir/graft/tensorflow"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Options struct {
	HiFiZMW      string
	AllZMW       string
	MinZMWDepth  int
	WindowRadius int
	ModelDir     string
	OutPrefix    string
}

func loadModel(modelPath string, modelNames []string) (*tf.SavedModel, error) {
	model, err := tf.LoadSavedModel(modelPath, modelNames, nil) // 载入模型
	if err != nil {
		log.Printf("LoadSavedModel err: %v", err)
		return nil, err
	}

	log.Println("list possible ops in graphs")
	for _, op := range model.Graph.Operations() {
		log.Printf("Op name: %v", op.Name())
	}

	return model, nil
}

func RunTask(opts Options) error {

	modelDir := opts.ModelDir
	hiFiZMW := opts.HiFiZMW
	allZMW := opts.AllZMW
	minzmwDepth := opts.MinZMWDepth
	windowRadius := opts.WindowRadius

	//step1: model loading
	model, err := loadModel(modelDir, []string{"serve"})
	if err != nil {
		log.Printf("loadModel err:%v, modelDir:%s", err, modelDir)
		return err
	}
	log.Printf("loaded model:%v", model)

	//step2: read in read-level prediction
	ZMWprobDictNorm, countBasedPredDict, depthOfMolCovDict := normZMWprediction_HCO(hiFiZMW, allZMW, minzmwDepth)
	XTest, XAveWeights, ZMWprobDictPosListSorted := CreateDataMatrixDistanceTesting(ZMWprobDictNorm, windowRadius)

	//step3: do the prediction
	probes, err := predict.Predict(model, XTest)
	if err != nil {
		log.Printf("model predict err:%+v", err)
		return err
	}
	//step4: format to write out
	outFileName := opts.OutPrefix + ".JointCpG.MeLoDe-Combo.txt"
	err = writeFile(outFileName, countBasedPredDict, depthOfMolCovDict, probes, XAveWeights, ZMWprobDictPosListSorted)
	if err != nil {
		log.Printf("writeFile err:%+v, outFileName:%s", err, outFileName)
		return err
	}
	return nil
}

func readFile(filename string) ([]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	gz, err := gzip.NewReader(f)
	if err != nil {
		return nil, err
	}
	defer gz.Close()

	var lines []string
	scanner := bufio.NewScanner(gz)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	return lines, scanner.Err()
}

func writeFile(filename string, countBasedPredDict map[string][]float32, depthOfMolCovDict map[string]int,
	probes []float32, XAveWeights []float64, ZMWprobDictPosListSorted []string) error {

	length := len(probes)
	if len(XAveWeights) != length || len(ZMWprobDictPosListSorted) != length {
		log.Printf("data length check failed, len(probes):%d,"+
			"len(XAveWeights):%d, len(ZMWprobDictPosListSorted):%d", len(probes),
			len(XAveWeights), len(ZMWprobDictPosListSorted))
		return fmt.Errorf("data length check failed")
	}
	file, err := os.Create(filename)
	if err != nil {
		log.Fatalf("could not open file %q:", err)
		return err
	}
	defer file.Close()

	// 创建一个写入器
	writer := bufio.NewWriter(file)

	var chrs, postions []string
	var depthOfMolCov []int
	var countBasedPred [][]float32
	for _, chrpos := range ZMWprobDictPosListSorted {
		chr, pos := parsePos(chrpos)
		chrs = append(chrs, chr)
		postions = append(postions, pos)
		depthOfMolCov = append(depthOfMolCov, depthOfMolCovDict[chrpos])
		countBasedPred = append(countBasedPred, countBasedPredDict[chrpos])
	}
	header := fmt.Sprintf("Chr\tPosC\tmolecule_cov\tmodel_based\tcount_based\tMe\tunMeth\tcisCpG_ave_weights\n")
	_, err = writer.WriteString(header)
	if err != nil {
		log.Println("Error writing line:", err)
		return err
	}
	for i := 0; i < length; i++ {
		line := fmt.Sprintf("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f",
			chrs[i], postions[i], depthOfMolCov[i], probes[i], countBasedPred[i][0], countBasedPred[i][1], countBasedPred[i][2], XAveWeights[i])
		_, err = writer.WriteString(line + "\n")
		if err != nil {
			log.Println("Error writing line:", err)
			return err
		}
	}
	// 将缓冲区的数据刷新到文件中
	err = writer.Flush()
	if err != nil {
		log.Println("Error flushing writer:", err)
		return err
	}
	return nil
}

func normZMWprediction_HCO(ZMWprefile_HiFi, ALLZMW_prefile string, minZMWDepth int) (map[string][]float32, map[string][]float32, map[string]int) {
	ZMWprobDictRaw := make(map[string][][]float32)
	ZMWprobDictNorm := make(map[string][]float32)
	countBasedPredDict := make(map[string][]float32)
	depthOfMolCovDict := make(map[string]int)

	processedPos := make(map[string]struct{})
	processedHiFiZMWID := make(map[string]struct{})

	// Read HiFi file
	hifiLines, err := readFile(ZMWprefile_HiFi)
	if err != nil {
		log.Fatal(err)
	}
	for _, line := range hifiLines {
		if !strings.Contains(line, "nan") && !strings.HasPrefix(line, "Chr") {
			fields := strings.Fields(line)
			Chr, PosC, ZMWid, ZMWprob := fields[0], fields[1], fields[2], fields[8]
			linekey := Chr + "_" + PosC
			zmwkey := linekey + "_" + ZMWid
			prob, _ := strconv.ParseFloat(ZMWprob, 64)
			if _, exists := ZMWprobDictRaw[linekey]; !exists {
				ZMWprobDictRaw[linekey] = [][]float32{{}, {}, {}}
			}
			ZMWprobDictRaw[linekey][0] = append(ZMWprobDictRaw[linekey][0], float32(prob))
			processedPos[linekey] = struct{}{}
			processedHiFiZMWID[zmwkey] = struct{}{}
		}
	}

	// Read closed and open file
	if ALLZMW_prefile != "NA" {
		closedOpenLines, err := readFile(ALLZMW_prefile)
		if err != nil {
			log.Fatal(err)
		}
		for _, line := range closedOpenLines {
			if !strings.Contains(line, "nan") && !strings.HasPrefix(line, "Chr") {
				fields := strings.Fields(line)
				Chr, PosC, ZMWid, ZMWtype, ZMWprob := fields[0], fields[1], fields[2], fields[3], fields[8]
				linekey := Chr + "_" + PosC
				zmwkey := linekey + "_" + ZMWid
				if _, exists := processedHiFiZMWID[zmwkey]; !exists {
					prob, _ := strconv.ParseFloat(ZMWprob, 64)
					if _, exists := ZMWprobDictRaw[linekey]; !exists {
						ZMWprobDictRaw[linekey] = [][]float32{{}, {}, {}}
					}
					if ZMWtype == "C" {
						ZMWprobDictRaw[linekey][1] = append(ZMWprobDictRaw[linekey][1], float32(prob))
					} else if ZMWtype == "O" {
						ZMWprobDictRaw[linekey][2] = append(ZMWprobDictRaw[linekey][2], float32(prob))
					}
					processedPos[linekey] = struct{}{}
				}
			}
		}
	}

	// Normalize probabilities
	for pos, lists := range ZMWprobDictRaw {
		if len(lists[0])+len(lists[1])+len(lists[2]) >= minZMWDepth {
			totalMolecules := len(lists[0]) + len(lists[1]) + len(lists[2])
			meMolecules := 0
			for _, prob := range lists[0] {
				if prob > 0.5 {
					meMolecules++
				}
			}
			for _, prob := range lists[1] {
				if prob > 0.5 {
					meMolecules++
				}
			}
			for _, prob := range lists[2] {
				if prob > 0.5 {
					meMolecules++
				}
			}
			countBasedDNAmFreq := float32(meMolecules) / float32(totalMolecules)
			countBasedPredDict[pos] = []float32{countBasedDNAmFreq, float32(meMolecules), float32(totalMolecules - meMolecules)}
			depthOfMolCovDict[pos] = totalMolecules

			HiFiNormProbList := calc.GetNormalizedHisto(lists[0])
			log.Printf("pos:%s, list[0]:%v, HiFiNormProbList:%+v", pos, lists[0], HiFiNormProbList)

			ClosedNormProbList := calc.GetNormalizedHisto(lists[1])
			OpenNormProbList := calc.GetNormalizedHisto(lists[2])

			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], HiFiNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], ClosedNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], OpenNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], countBasedDNAmFreq)
			log.Printf("ZMWprobDictNorm[%s]:%+v", pos, ZMWprobDictNorm[pos])
			log.Printf("countBasedPredDict[%s]:%+v", pos, countBasedPredDict[pos])
		}
	}

	return ZMWprobDictNorm, countBasedPredDict, depthOfMolCovDict
}

func CreateDataMatrixDistanceTesting(ZMWprobDictNorm map[string][]float32, R int) ([][][]float32, []float64, []string) {
	totalNormHistos := [][]float32{}
	totalPos := []int{}
	totalAveWeights := []float64{}
	ZMWprobDictPosList := []string{}

	for pos := range ZMWprobDictNorm {
		ZMWprobDictPosList = append(ZMWprobDictPosList, pos)
	}

	sort.Slice(ZMWprobDictPosList, func(i, j int) bool {
		iPos, _ := strconv.Atoi(strings.Split(ZMWprobDictPosList[i], "_")[1])
		jPos, _ := strconv.Atoi(strings.Split(ZMWprobDictPosList[j], "_")[1])
		return iPos < jPos
	})

	for _, pos := range ZMWprobDictPosList {
		totalNormHistos = append(totalNormHistos, ZMWprobDictNorm[pos])
		_, posC := parsePos(pos)
		intPosC, _ := strconv.Atoi(posC)
		totalPos = append(totalPos, intPosC)
	}

	featPad := totalNormHistos
	featuresWindow := calc.SlidingWindowView2DArray(featPad, 2*R+1, 0)
	N := len(featuresWindow)

	featDis := totalPos
	distanceWin := calc.SlidingWindowViewArray(featDis, 2*R+1)

	centeredDistances := [][]float64{}
	for i := 0; i < N; i++ {
		window := distanceWin[i]
		centeredValue := window[len(window)/2]
		distances := make([]float64, len(window))
		for j, val := range window {
			distances[j] = math.Abs(float64(val - centeredValue))
		}

		distancesV1 := make([]float64, len(distances))
		distancesV2 := make([]float64, len(distances))
		for j, val := range distances {
			distancesV1[j] = 100.0 / (math.Sqrt(val) + 100.0)
			distancesV2[j] = (distancesV1[j] - 0.5) / 0.5
		}
		centeredDistances = append(centeredDistances, distancesV2)
	}

	M := len(featuresWindow[0][0]) //M = 31
	featuresWindowMerge := make([][][]float32, N)
	for i := range featuresWindowMerge {
		featuresWindowMerge[i] = make([][]float32, 2*R+1)
		for j := range featuresWindowMerge[i] {
			featuresWindowMerge[i][j] = make([]float32, M+1) //N x (2R+1) x (M+1)
		}
	}

	for i := 0; i < N; i++ {
		a := centeredDistances[i]
		aveA := calc.Mean(a)
		totalAveWeights = append(totalAveWeights, aveA)
		b := featuresWindow[i]

		for j := 0; j < 2*R+1; j++ {
			copy(featuresWindowMerge[i][j], b[j])
			featuresWindowMerge[i][j][M] = float32(a[j])
		}
	}

	return featuresWindowMerge, totalAveWeights, ZMWprobDictPosList[R : len(ZMWprobDictPosList)-R]
}

// Helper function to parse position string
func parsePos(pos string) (string, string) {
	parts := strings.Split(pos, "_")
	return parts[0], parts[1]
}
