package task

import (
	"bufio"
	"fmt"
	"github.com/Taichidasheen/tag_decode/pkg/record_flag"
	"github.com/Taichidasheen/tag_decode/pkg/record_tag"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog/log"
	"io"
	"os"
	"strconv"
	"strings"
	"sync"
)

func RunDecodeTask(opts Options) error {

	bamFile := opts.BamFile
	outFile := opts.OutFile
	topN := opts.TopN

	var wg sync.WaitGroup

	wg.Add(3)
	recordChan := make(chan *sam.Record, 1000)
	resultChan := make(chan string, 3000)
	//predictResultChan := make(chan *sam.Record, 1000)

	//读取bam文件
	go readBam(&wg, recordChan, bamFile, topN)

	//异步处理record
	go processAlignedHiFiRecord(&wg, recordChan, resultChan)

	//resultPath := "/storage/yangjianLab/westlakechat/subreads_locate/result.txt"
	//异步写入result
	go writeTextResult(&wg, outFile, resultChan)

	wg.Wait()

	log.Info().Msgf("exit...")

	return nil

}

func readBam(wg *sync.WaitGroup, recordChan chan *sam.Record, bamFilePath string, topN int) {
	defer wg.Done()
	defer close(recordChan)

	//bamFilePath := "/storage/yangjianLab/caoyujie/project/meth/Bam_file/4_hifi_subreads/HG01109_WT_hifi/5x.sort.HG01109_WT_hifi.aln.bam"
	bamFile, err := os.Open(bamFilePath)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	defer bamFile.Close()
	ok, err := bgzf.HasEOF(bamFile)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	if !ok {
		log.Printf("file has no bgzf magic block: may be truncated")
		return
	}
	//读写并发度
	//concurrency := 8
	concurrency := 5
	//bam reader
	bamReader, err := bam.NewReader(bamFile, concurrency)
	if err != nil {
		log.Error().Msgf("could not read bam:%v", err)
		return
	}
	defer bamReader.Close()

	var count int
	fmt.Printf("begin count...")
	for {
		record, err := bamReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Error().Msgf("error reading bam: %v", err)
			return
		}
		count++
		if count%10000 == 0 {
			log.Printf("count:%d", count)
		}
		if topN != 0 && count == topN+1 {
			log.Printf("has processed topN:%d rows, break...", topN)
			break
		}
		recordChan <- record
	}
}

func processAlignedHiFiRecord(wg *sync.WaitGroup, recordChan chan *sam.Record,
	resultChan chan string) {
	defer wg.Done()

	//记得close resultChan, 否则会deadlock
	defer close(resultChan)

	for record := range recordChan {

		log.Debug().Msgf("processing readname:%s", record.Name)

		if !record_flag.IsSecondary(record.Flags) && !record_flag.IsSupplementary(record.Flags) {
			//mm
			mmTag, err := record_tag.ExtractMMTag(record)
			if err != nil {
				log.Error().Msgf("extract mm tag err, record name:%s", record.Name)
				continue
			}
			//ml
			mlAux, exist := record.Tag([]byte{'M', 'L'})
			if !exist {
				log.Error().Msgf("record ml not exist, record name:%s", record.Name)
				continue
			}
			mlTag, ok := mlAux.Value().([]uint8)
			if !ok {
				log.Error().Msgf("record ml invalid, record name:%s", record.Name)
				continue
			}
			readSeqList := record.Seq.Expand()
			alnRefChr := record.Ref.Name()
			alnRefStart := record.Pos
			readIsReverse := record_flag.IsReverse(record.Flags)
			readCigar := record.Cigar
			readName := record.Name
			parts := strings.Split(readName, "/")
			zmwName := parts[0] + "_" + parts[1]
			zmwType := "H"

			recordTag, err := record_tag.ExtractRecordTag(record)
			if err != nil {
				log.Error().Msgf("record_tag.ExtractRecordTag err:%v, record name:%s", err, record.Name)
				continue
			}
			//检查mm和ml的长度是否相等
			if len(mmTag) != len(mlTag) {
				log.Error().Msgf("mm and ml length not equal err:%v, record name:%s", err, record.Name)
				continue
			}

			//fn, rn
			fn := recordTag.Fn
			rn := recordTag.Rn
			totalSubreadsDep := "X"
			if fn != "X" && rn != "X" {
				fnInt, _ := strconv.Atoi(fn)
				rnInt, _ := strconv.Atoi(rn)
				totalSubreadsDep = fmt.Sprintf("%d", fnInt+rnInt)
			}
			haplotype := recordTag.HP
			haplotypeBlock := recordTag.PS
			strandSign := "+"

			tagOnSEQPos, tagOnSEQProb := mmMLTagOnSEQ(string(readSeqList), readIsReverse, mmTag, mlTag)
			tagOnRefPos, tagOnSEQProbUpdate := ConPosSEQToPosRefQuick(tagOnSEQPos, alnRefStart, readCigar, tagOnSEQProb)
			for i, refPos := range tagOnRefPos {
				var bi int
				if tagOnSEQProbUpdate[i] >= 0.5 {
					bi = 1
				} else {
					bi = 0
				}
				line := fmt.Sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%v\t%d",
					alnRefChr, refPos, zmwName, zmwType, totalSubreadsDep, haplotype, haplotypeBlock, strandSign, tagOnSEQProbUpdate[i], bi)
				resultChan <- line
			}
		}
	}

}

// ConPosSEQToPosRef 将 CpG 位点从读取序列上的位置映射到参考基因组的位置
func ConPosSEQToPosRef(tagOnSEQPos []int, alnRefStart int, readCigar sam.Cigar, tagOnSEQProb []float64) ([]int, []float64) {
	// 初始化输出
	tagOnRefPos := make([]string, len(tagOnSEQPos)) // 用字符串表示未映射的点为 "NA"
	for i := range tagOnRefPos {
		tagOnRefPos[i] = "NA"
	}

	// 遍历每个 CpG 位点
	for index, seqPos := range tagOnSEQPos {
		readPosBlkStart, readPosBlkEnd := 0, 0
		refWalkingPosStart := alnRefStart
		refWalkingPosEnd := 0
		flag := 1

		// 遍历 CIGAR 块
		for _, cigar := range readCigar {
			op := cigar.Type()
			count := cigar.Len()

			// 更新读取序列的范围
			if op == 0 || op == 1 || op == 4 || op == 7 || op == 8 { // M, I, S, =, X
				readPosBlkEnd += count
				readPosBlkStart = readPosBlkEnd - count
			}

			// 更新参考基因组的范围
			if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 { // M, D, N, =, X
				if refWalkingPosStart == alnRefStart && flag == 1 {
					refWalkingPosEnd = refWalkingPosStart + count
					flag = 2
				} else {
					refWalkingPosStart = refWalkingPosEnd
					refWalkingPosEnd += count
				}
			}

			// 检查 CpG 位点是否在读取序列范围内
			if readPosBlkStart <= seqPos && seqPos <= readPosBlkEnd {
				if op == 0 || op == 1 || op == 4 || op == 7 || op == 8 { // M, I, S, =, X
					lastOpNeeded := seqPos - readPosBlkStart
					if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 { // M, D, N, =, X
						refPos := refWalkingPosStart + lastOpNeeded
						tagOnRefPos[index] = fmt.Sprintf("%d", refPos)
						break
					}
				}
			}
		}
	}

	// 移除未映射的位点
	var tagOnRefPosUpdate []int
	var tagOnSEQProbUpdate []float64
	for i, pos := range tagOnRefPos {
		if pos != "NA" {
			refPos := 0
			fmt.Sscanf(pos, "%d", &refPos)
			tagOnRefPosUpdate = append(tagOnRefPosUpdate, refPos)
			tagOnSEQProbUpdate = append(tagOnSEQProbUpdate, tagOnSEQProb[i])
		}
	}

	return tagOnRefPosUpdate, tagOnSEQProbUpdate
}

// ConPosSEQToPosRefQuick 将 CpG 位点从读取序列上的位置映射到参考基因组的位置
func ConPosSEQToPosRefQuick(tagOnSEQPos []int, alnRefStart int, readCigar sam.Cigar, tagOnSEQProb []float64) ([]int, []float64) {
	// 初始化输出
	tagOnRefPos := make([]string, len(tagOnSEQPos)) // 用字符串表示未映射的点为 "NA"
	for i := range tagOnRefPos {
		tagOnRefPos[i] = "NA"
	}

	cpgBeginIdx := 0
	readPosBlkStart, readPosBlkEnd := 0, 0
	refWalkingPosStart := alnRefStart
	refWalkingPosEnd := 0
	flag := 1

	// 遍历 CIGAR 块
	for _, cigar := range readCigar {
		op := cigar.Type()
		count := cigar.Len()

		// 更新读取序列的范围
		if op == 0 || op == 1 || op == 4 || op == 7 || op == 8 { // M, I, S, =, X
			readPosBlkEnd += count
			readPosBlkStart = readPosBlkEnd - count
		}

		// 更新参考基因组的范围
		if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 { // M, D, N, =, X
			if refWalkingPosStart == alnRefStart && flag == 1 {
				refWalkingPosEnd = refWalkingPosStart + count
				flag = 2
			} else {
				refWalkingPosStart = refWalkingPosEnd
				refWalkingPosEnd += count
			}
		}

		//查看readPosBlkStart到readPosBlkEnd范围内有哪些read CpG位点
		//只有当op为M,=,X的时候才去执行read和ref的比对逻辑
		//也就是说只有当read和ref上的位置都向前移动的时候，才进行read和ref的比对
		if op == 0 || op == 7 || op == 8 {
			matchedCpgs, matchedCpgIdxs, nextIdx := findCpgPosOnReadInRange(tagOnSEQPos, cpgBeginIdx, readPosBlkStart, readPosBlkEnd)
			cpgBeginIdx = nextIdx
			for i, cpg := range matchedCpgs {
				lastOpNeeded := cpg - readPosBlkStart
				//refPos := refWalkingPosStart + lastOpNeeded
				refPos := refWalkingPosStart + lastOpNeeded + 1 //+1之后和python版本对上
				index := matchedCpgIdxs[i]
				tagOnRefPos[index] = fmt.Sprintf("%d", refPos)
			}
			if nextIdx >= len(tagOnSEQPos) {
				//对比结束
				break
			}
		}
	}

	// 移除未映射的位点
	var tagOnRefPosUpdate []int
	var tagOnSEQProbUpdate []float64
	for i, pos := range tagOnRefPos {
		if pos != "NA" {
			refPos := 0
			fmt.Sscanf(pos, "%d", &refPos)
			tagOnRefPosUpdate = append(tagOnRefPosUpdate, refPos)
			tagOnSEQProbUpdate = append(tagOnSEQProbUpdate, tagOnSEQProb[i])
		}
	}

	return tagOnRefPosUpdate, tagOnSEQProbUpdate
}

// 判断哪些cpg在当前给定的范围内，返回这些cpg的值
// cpgBeginIdx表示从tagOnSEQPos的哪个位置开始比对（cpgBeginIdx之前的位置可以不用比对）
func findCpgPosOnReadInRange(tagOnSEQPos []int, cpgBeginIdx int, readPosBlkStart, readPosBlkEnd int) ([]int, []int, int) {
	var matchedCpgs []int
	var matchedCpgIdxs []int
	var cpgNextIndex int //下一次查找overlappingCpg数组的开始位置

	for i := cpgBeginIdx; i < len(tagOnSEQPos); i++ {
		if readPosBlkStart <= tagOnSEQPos[i] && tagOnSEQPos[i] <= readPosBlkEnd {
			matchedCpgs = append(matchedCpgs, tagOnSEQPos[i])
			matchedCpgIdxs = append(matchedCpgIdxs, i)
			cpgNextIndex = i + 1
		}
		if readPosBlkStart > tagOnSEQPos[i] {
			cpgNextIndex = i + 1
		}
		if readPosBlkEnd < tagOnSEQPos[i] {
			//不用再比对了，等到refstart和refend更新后再从i位置开始比对
			cpgNextIndex = i
			break
		}
	}
	return matchedCpgs, matchedCpgIdxs, cpgNextIndex
}

// mmMLTagOnSEQ 计算 CpG 位点位置和甲基化概率
func mmMLTagOnSEQ(readSEQString string, readIsReverse bool, mmtag []int, mltag []uint8) ([]int, []float64) {
	readQueryLength := len(readSEQString)

	var tagOnTOPPos []int
	var tagOnTOPProb []float64

	var tagOnSEQPos []int
	var tagOnSEQProb []float64

	// 获取正义链
	var topStrand string
	if readIsReverse {
		reverseSeq := reverseString(readSEQString)
		topStrand = string(seqComplementary([]byte(reverseSeq)))
	} else {
		topStrand = readSEQString
	}

	// 获取正义链上的C位点位置 (1-based)
	var CposOnTopStrand []int
	for i, base := range topStrand {
		if base == 'C' {
			CposOnTopStrand = append(CposOnTopStrand, i+1)
		}
	}

	log.Debug().Msgf("mmMLTagOnSEQ mmtag:%+v,len(mmtag):%d", mmtag, len(mmtag))
	log.Debug().Msgf("mmMLTagOnSEQ mltag:%+v,len(mltag):%d", mltag, len(mltag))
	cIndex := 0
	for i := 0; i < len(mmtag); i++ {
		cIndex += mmtag[i] + 1                    //当前cg位点中的c是第几个c。c的计数从1开始
		cpgPos := CposOnTopStrand[cIndex-1]       //CposOnTopStrand，1-based
		tagOnTOPPos = append(tagOnTOPPos, cpgPos) //tagOnTOPPos，1-based
		tagOnTOPProb = append(tagOnTOPProb, float64(mltag[i])/256)
	}

	// 如果是反向链，转换位点位置并翻转概率
	if readIsReverse {
		for i := len(tagOnTOPPos) - 1; i >= 0; i-- {
			tagOnSEQPos = append(tagOnSEQPos, readQueryLength-tagOnTOPPos[i]-1)
			tagOnSEQProb = append(tagOnSEQProb, tagOnTOPProb[i])
		}
	} else {
		for _, pos := range tagOnTOPPos {
			tagOnSEQPos = append(tagOnSEQPos, pos-1) //tagOnSEQPos， 0-based
		}
		tagOnSEQProb = tagOnTOPProb
	}

	return tagOnSEQPos, tagOnSEQProb
}

// 计算反向互补序列
func seqComplementary(seq []byte) []byte {
	complementedSeq := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		switch seq[i] {
		case 'A':
			complementedSeq[i] = 'T'
		case 'T':
			complementedSeq[i] = 'A'
		case 'C':
			complementedSeq[i] = 'G'
		case 'G':
			complementedSeq[i] = 'C'
		default:
			log.Warn().Msgf("invalid bytes:%v", seq[i])
		}
	}
	return complementedSeq
}

// reverseString 翻转字符串
func reverseString(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func writeTextResult(wg *sync.WaitGroup, outFile string, resultChan chan string) {
	defer wg.Done()

	file, err := os.Create(outFile)
	if err != nil {
		log.Error().Msgf("could not open file %q:", err)
		return
	}
	defer file.Close()

	// 创建一个写入器
	writer := bufio.NewWriter(file)
	for line := range resultChan {
		_, err := writer.WriteString(line + "\n")
		if err != nil {
			log.Error().Msgf("Error writing line:%v", err)
			return
		}
	}
	// 将缓冲区的数据刷新到文件中
	err = writer.Flush()
	if err != nil {
		log.Error().Msgf("Error flushing writer:%v", err)
		return
	}
}
