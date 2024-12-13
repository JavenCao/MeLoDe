package main

import (
	"flag"
	"fmt"
	"github.com/Taichidasheen/site_predict/pkg/task"
	"log"
	"os"
	"runtime"
	"runtime/pprof"
)

func main() {

	//input param
	var hifiZMW string
	var minZMWDepth int
	var windowRadius int
	var modelDir string
	var outPrefix string

	//processing param
	var processor int

	//debug param
	var cpuprofile string
	var topN int

	//input related parameters
	flag.StringVar(&hifiZMW, "HiFiZMWpre", "", "HiFi reads - single molecule level prediction. File format: txt.gz")
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

	flag.Parse()

	fmt.Println("HiFiZMWpre:", hifiZMW)
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
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	log.Printf("begin to process site predict task")

	opts := task.Options{
		HiFiZMW:      hifiZMW,
		MinZMWDepth:  minZMWDepth,
		WindowRadius: windowRadius,
		ModelDir:     modelDir,
		OutPrefix:    outPrefix,
	}

	err := task.RunHiFiTask(opts)
	if err != nil {
		log.Fatalf("run task err:%+v, opts:%+v", err, opts)
		return
	}

	log.Printf("site predict task is completed!")

}
