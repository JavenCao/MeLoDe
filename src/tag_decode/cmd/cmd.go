package main

import (
	"flag"
	"github.com/Taichidasheen/tag_decode/pkg/task"
	"github.com/rs/zerolog"
	"github.com/rs/zerolog/log"
	"os"
	"runtime"
	"runtime/pprof"
)

func main() {
	var bamFilePath string
	var outFile string

	//debug param
	var cpuprofile string
	var topN int
	var logLevel int

	flag.StringVar(&bamFilePath, "bamfile", "", "Aligned HiFi Bam file with MM/ML tags")
	flag.StringVar(&outFile, "outfile", "", "Read levle prediction file")

	//debug parameters
	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to this file")
	flag.IntVar(&topN, "topN", 0, "just process top N rows")
	flag.IntVar(&logLevel, "loglevel", 1, "0-debug,1-info,2-warn,3-error")

	flag.Parse()
	log.Output(os.Stdout)
	zerolog.SetGlobalLevel(zerolog.Level(logLevel))

	log.Info().Msgf("cpuprofile:%v", cpuprofile)
	log.Info().Msgf("topN:%v", topN)
	log.Info().Msgf("loglevel:%v", logLevel)

	log.Info().Msgf("maxProcs:%v", runtime.GOMAXPROCS(0))

	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Err(err)
			return
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	opts := task.Options{
		BamFile: bamFilePath,
		OutFile: outFile,
		TopN:    topN,
	}

	err := task.RunDecodeTask(opts)
	if err != nil {
		log.Error().Msgf("run task err:%+v, opts:%+v", err, opts)
		return
	}

	log.Info().Msgf("tag decode task is completed!")

}
