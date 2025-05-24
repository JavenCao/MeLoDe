package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strings"
)

type OrderedMap struct {
	keys []string
	data map[string]string
}

func NewOrderedMap() *OrderedMap {
	return &OrderedMap{
		keys: []string{},
		data: make(map[string]string),
	}
}

func (om *OrderedMap) Set(key, value string) {
	if _, exists := om.data[key]; !exists {
		om.keys = append(om.keys, key)
	}
	om.data[key] = value
}

func (om *OrderedMap) Has(key string) bool {
	_, exists := om.data[key]
	return exists
}

func (om *OrderedMap) MergeFrom(other *OrderedMap) {
	for _, key := range other.keys {
		if !om.Has(key) {
			om.Set(key, other.data[key])
		}
	}
}

func (om *OrderedMap) WriteToFile(filename string) error {
	f, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer f.Close()

	for _, key := range om.keys {
		_, err := f.WriteString(om.data[key] + "\n")
		if err != nil {
			return err
		}
	}
	return nil
}

func readFileToOrderedMap(filename string, exclude *OrderedMap) (*OrderedMap, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	omap := NewOrderedMap()
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) < 2 {
			continue
		}
		key := fields[0] + "_" + fields[1]
		if exclude == nil || !exclude.Has(key) {
			omap.Set(key, line)
		}
	}
	return omap, nil
}

func main() {
	// Parse args
	Hfile := flag.String("Hfile", "", "site-level file (H_model)")
	HCfile := flag.String("HCfile", "", "site-level file (HC_model)")
	HCOfile := flag.String("HCOfile", "", "site-level file (HCO_model)")
	output := flag.String("o", "", "output file prefix")
	flag.Parse()

	if *Hfile == "" || *HCfile == "" || *HCOfile == "" || *output == "" {
		flag.Usage()
		os.Exit(1)
	}

	outH := *output + ".MeLoDe_H.site.txt"
	outHHC := *output + ".MeLoDe_H_HC.txt"
	outHHCHCO := *output + ".MeLoDe_H_HC_HCO.txt"

	// Step 1: Read H
	HMap, err := readFileToOrderedMap(*Hfile, nil)
	if err != nil {
		fmt.Println("Error reading H file:", err)
		return
	}

	// Step 2: Read HC (exclude keys in H)
	HCMap, err := readFileToOrderedMap(*HCfile, HMap)
	if err != nil {
		fmt.Println("Error reading HC file:", err)
		return
	}

	// Step 3: Read HCO (exclude keys in H and HC)
	excludeHHC := NewOrderedMap()
	excludeHHC.MergeFrom(HMap)
	excludeHHC.MergeFrom(HCMap)

	HCOMap, err := readFileToOrderedMap(*HCOfile, excludeHHC)
	if err != nil {
		fmt.Println("Error reading HCO file:", err)
		return
	}

	// Step 4: Merge
	HHCMap := NewOrderedMap()
	HHCMap.MergeFrom(HMap)
	HHCMap.MergeFrom(HCMap)

	HHCHCOMap := NewOrderedMap()
	HHCHCOMap.MergeFrom(HHCMap)
	HHCHCOMap.MergeFrom(HCOMap)

	// Step 5: Write output files
	if err := HMap.WriteToFile(outH); err != nil {
		fmt.Println("Error writing H file:", err)
	}
	if err := HHCMap.WriteToFile(outHHC); err != nil {
		fmt.Println("Error writing H_HC file:", err)
	}
	if err := HHCHCOMap.WriteToFile(outHHCHCO); err != nil {
		fmt.Println("Error writing H_HC_HCO file:", err)
	}
}
