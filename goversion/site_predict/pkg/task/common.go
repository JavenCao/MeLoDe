package task

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"github.com/Taichidasheen/site_predict/pkg/calc"
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
