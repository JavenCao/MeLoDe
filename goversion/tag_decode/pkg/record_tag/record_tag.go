package record_tag

import (
	"fmt"
	"github.com/biogo/hts/sam"
	"github.com/rs/zerolog/log"
	"strconv"
	"strings"
)

// RecordTag Notice:fn,rn这些tag在bam文件中是fn,rn(区分大小写)
type RecordTag struct {
	Fn string
	Rn string
	// HP:i:2   PS:i:11457
	HP string
	PS string
}

func ExtractRecordTag(record *sam.Record) (*RecordTag, error) {

	var fn, rn string
	//fn
	fnAux, exist := record.Tag([]byte{'f', 'n'})
	if !exist {
		fn = "X"
	} else {
		fnVal, ok := fnAux.Value().(int32)
		if !ok {
			log.Error().Msgf("record fn invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("fn invalid")
		}
		fn = fmt.Sprintf("%d", fnVal)
	}

	//rn
	rnAux, exist := record.Tag([]byte{'r', 'n'})
	if !exist {
		rn = "X"
	} else {
		rnVal, ok := rnAux.Value().(int32)
		if !ok {
			log.Error().Msgf("record rn invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("rn invalid")
		}
		rn = fmt.Sprintf("%d", rnVal)
	}

	//fmt.Println("record name:", record.Name, "fn:", fn, " rn:", rn)

	// HP:i:2   PS:i:11457
	var HP, PS string
	HPAux, exist := record.Tag([]byte{'H', 'P'})
	if !exist {
		//log.Error().Msgf("record HP not exist, record name:%s", record.Name)
		HP = "X"
	} else {
		HPVal, ok := HPAux.Value().(int32)
		if !ok {
			log.Error().Msgf("record HP invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("HP invalid")
		}
		HP = fmt.Sprintf("%d", HPVal)
	}

	PSAux, exist := record.Tag([]byte{'P', 'S'})
	if !exist {
		//log.Error().Msgf("record PS not exist, record name:%s", record.Name)
		PS = "X"
	} else {
		PSVal, ok := PSAux.Value().(int32)
		if !ok {
			log.Error().Msgf("record PS invalid, record name:%s", record.Name)
			return nil, fmt.Errorf("PS invalid")
		}
		PS = fmt.Sprintf("%d", PSVal)
	}

	recordTag := &RecordTag{
		Fn: fn,
		Rn: rn,
		HP: HP,
		PS: PS,
	}
	return recordTag, nil
}

func ExtractMMTag(record *sam.Record) ([]int, error) {
	mmAux, exist := record.Tag([]byte{'M', 'M'})
	if !exist {
		log.Error().Msgf("record mm not exist, record name:%s", record.Name)
		return nil, fmt.Errorf("mm tag not exist")
	}
	mmTagStr := mmAux.String()
	log.Debug().Msgf("mmTagStr:%+v", mmTagStr)
	if !strings.HasSuffix(mmTagStr, ";") {
		log.Error().Msgf("record mm invalid, record name:%s", record.Name)
		return nil, fmt.Errorf("mm tag invalid")
	}
	mmTagStr = mmTagStr[0 : len(mmTagStr)-1] //去掉末尾的;
	mmTagStrArr := strings.Split(mmTagStr, ",")
	var mmTagArr []int
	for i := 1; i < len(mmTagStrArr); i++ {
		mmstr := mmTagStrArr[i]
		mmItem, err := strconv.Atoi(mmstr)
		if err != nil {
			log.Error().Msgf("record mm strconv.Atoi err:%+v, record name:%s", err, record.Name)
			return nil, err
		}
		mmTagArr = append(mmTagArr, mmItem)
	}
	return mmTagArr, nil
}
