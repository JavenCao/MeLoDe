package task

import (
	"github.com/Taichidasheen/site_predict/pkg/calc"
	"github.com/Taichidasheen/site_predict/pkg/predict"
	"github.com/rs/zerolog/log"
	"strconv"
	"strings"
)

func RunComboTaskHCO(opts Options) error {

	modelDir := opts.ModelDir
	hiFiZMW := opts.HiFiZMW
	allZMW := opts.AllZMW
	minzmwDepth := opts.MinZMWDepth
	windowRadius := opts.WindowRadius

	//step1: model loading
	model, err := loadModel(modelDir, []string{"serve"})
	if err != nil {
		log.Error().Msgf("loadModel err:%v, modelDir:%s", err, modelDir)
		return err
	}
	log.Info().Msgf("loaded model:%v", model)

	//step2: read in read-level prediction
	ZMWprobDictNorm, countBasedPredDict, depthOfMolCovDict := normZMWprediction_HCO(hiFiZMW, allZMW, minzmwDepth)
	XTest, XAveWeights, ZMWprobDictPosListSorted := CreateDataMatrixDistanceTesting(ZMWprobDictNorm, windowRadius)

	//step3: do the prediction
	probes, err := predict.Predict(model, XTest)
	if err != nil {
		log.Error().Msgf("model predict err:%+v", err)
		return err
	}
	//step4: format to write out
	outFileName := opts.OutPrefix + ".JointCpG.MeLoDe-Combo.HCO.txt"
	err = writeFile(outFileName, countBasedPredDict, depthOfMolCovDict, probes, XAveWeights, ZMWprobDictPosListSorted)
	if err != nil {
		log.Error().Msgf("writeFile err:%+v, outFileName:%s", err, outFileName)
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
		log.Fatal().Msgf("%s readFile err:%+v", ZMWprefile_HiFi, err)
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
			log.Fatal().Msgf("%s readFile err:%+v", ALLZMW_prefile, err)
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
			log.Debug().Msgf("pos:%s, list[0]:%v, HiFiNormProbList:%+v", pos, lists[0], HiFiNormProbList)

			ClosedNormProbList := calc.GetNormalizedHisto(lists[1])
			OpenNormProbList := calc.GetNormalizedHisto(lists[2])

			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], HiFiNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], ClosedNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], OpenNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], countBasedDNAmFreq)
			log.Debug().Msgf("ZMWprobDictNorm[%s]:%+v", pos, ZMWprobDictNorm[pos])
			log.Debug().Msgf("countBasedPredDict[%s]:%+v", pos, countBasedPredDict[pos])
		}
	}

	return ZMWprobDictNorm, countBasedPredDict, depthOfMolCovDict
}

func RunComboTaskHC(opts Options) error {

	modelDir := opts.ModelDir
	hiFiZMW := opts.HiFiZMW
	allZMW := opts.AllZMW
	minzmwDepth := opts.MinZMWDepth
	windowRadius := opts.WindowRadius

	//step1: model loading
	model, err := loadModel(modelDir, []string{"serve"})
	if err != nil {
		log.Error().Msgf("loadModel err:%v, modelDir:%s", err, modelDir)
		return err
	}
	log.Info().Msgf("loaded model:%v", model)

	//step2: read in read-level prediction
	ZMWprobDictNorm, countBasedPredDict, depthOfMolCovDict := normZMWprediction_HC(hiFiZMW, allZMW, minzmwDepth)
	XTest, XAveWeights, ZMWprobDictPosListSorted := CreateDataMatrixDistanceTesting(ZMWprobDictNorm, windowRadius)

	//step3: do the prediction
	probes, err := predict.Predict(model, XTest)
	if err != nil {
		log.Error().Msgf("model predict err:%+v", err)
		return err
	}
	//step4: format to write out
	outFileName := opts.OutPrefix + ".JointCpG.MeLoDe-Combo.HC.txt"
	err = writeFile(outFileName, countBasedPredDict, depthOfMolCovDict, probes, XAveWeights, ZMWprobDictPosListSorted)
	if err != nil {
		log.Error().Msgf("writeFile err:%+v, outFileName:%s", err, outFileName)
		return err
	}
	return nil
}

func normZMWprediction_HC(ZMWprefile_HiFi, ALLZMW_prefile string, minZMWDepth int) (map[string][]float32, map[string][]float32, map[string]int) {
	ZMWprobDictRaw := make(map[string][][]float32)
	ZMWprobDictNorm := make(map[string][]float32)
	countBasedPredDict := make(map[string][]float32)
	depthOfMolCovDict := make(map[string]int)

	processedPos := make(map[string]struct{})
	processedHiFiZMWID := make(map[string]struct{})

	// Read HiFi file
	hifiLines, err := readFile(ZMWprefile_HiFi)
	if err != nil {
		log.Fatal().Msgf("%s readFile err:%+v", ZMWprefile_HiFi, err)
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
			log.Fatal().Msgf("%s readFile err:%+v", ALLZMW_prefile, err)
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
						//skip type O
						continue
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
			log.Debug().Msgf("pos:%s, list[0]:%v, HiFiNormProbList:%+v", pos, lists[0], HiFiNormProbList)

			ClosedNormProbList := calc.GetNormalizedHisto(lists[1])
			OpenNormProbList := calc.GetNormalizedHisto(lists[2])

			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], HiFiNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], ClosedNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], OpenNormProbList...)
			ZMWprobDictNorm[pos] = append(ZMWprobDictNorm[pos], countBasedDNAmFreq)
			log.Debug().Msgf("ZMWprobDictNorm[%s]:%+v", pos, ZMWprobDictNorm[pos])
			log.Debug().Msgf("countBasedPredDict[%s]:%+v", pos, countBasedPredDict[pos])
		}
	}

	return ZMWprobDictNorm, countBasedPredDict, depthOfMolCovDict
}
