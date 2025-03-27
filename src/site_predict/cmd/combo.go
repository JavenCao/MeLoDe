package main

import (
	"flag"
	"fmt"
	"github.com/Taichidasheen/site_predict/pkg/task"
	"github.com/rs/zerolog"
	"github.com/rs/zerolog/log"
	"os"
	"runtime"
	"runtime/pprof"
)

func main() {

	//input param
	var hifiZMW string
	var subreadsZMW string
	var minZMWDepth int
	var windowRadius int
	var modelDir string
	var outPrefix string

	//processing param
	var processor int

	//debug param
	var cpuprofile string
	var topN int
	var logLevel int

	//input related parameters
	flag.StringVar(&hifiZMW, "HiFiZMWpre", "", "HiFi reads - single molecule level prediction. File format: txt.gz")
	flag.StringVar(&subreadsZMW, "SubreadsZMWpre", "NA", "ZMW prediction using subreads. File format: txt.gz")
	flag.IntVar(&minZMWDepth, "minZMWdepth", 4, "min ZMW_molecule/HiFi_reads/CCS_reads depth")
	flag.IntVar(&windowRadius, "windowRadius", 5, "how many CpG sites in a sliding window? (2 * r + 1)")
	flag.StringVar(&modelDir, "modelDir", "", "Joint-CpG model directory for MeLoDe-Combo")

	// processing parameters
	flag.IntVar(&processor, "p", 0, "Parallelism processors")

	// output related parameters
	flag.StringVar(&outPrefix, "outputfile", "", "CpG site level methylation frequency")

	//debug parameters
	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to this file")
	flag.IntVar(&topN, "topN", 0, "just process top N rows")
	flag.IntVar(&logLevel, "loglevel", 3, "0-debug,1-info,2-warn,3-error")

	flag.Parse()
	log.Output(os.Stdout)
	zerolog.SetGlobalLevel(zerolog.Level(logLevel))

	fmt.Println("HiFiZMWpre:", hifiZMW)
	fmt.Println("SubreadsZMWpre:", subreadsZMW)
	fmt.Println("minZMWdepth:", minZMWDepth)
	fmt.Println("windowRadius:", windowRadius)
	fmt.Println("modelDir:", modelDir)
	fmt.Println("o:", outPrefix)

	fmt.Println("processor:", processor)

	fmt.Println("cpuprofile:", cpuprofile)
	fmt.Println("topN:", topN)

	fmt.Println("maxProcs:", runtime.GOMAXPROCS(0))

	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Fatal().Msgf("err:%v", err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	log.Printf("begin to process site predict task")

	opts := task.Options{
		HiFiZMW:      hifiZMW,
		AllZMW:       subreadsZMW,
		MinZMWDepth:  minZMWDepth,
		WindowRadius: windowRadius,
		ModelDir:     modelDir,
		OutPrefix:    outPrefix,
	}

	err := task.RunComboTask(opts)
	if err != nil {
		log.Fatal().Msgf("run task err:%+v, opts:%+v", err, opts)
		return
	}

	log.Printf("site predict task is completed!")

}
