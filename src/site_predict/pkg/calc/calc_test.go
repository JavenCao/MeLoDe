package calc

import (
	"testing"
)

func TestGetNormalizedHisto(t *testing.T) {

	//probs := []float64{0.1, 0.2, 0.3, 0.4, 0.5}
	//t.Logf("Histogram:%+v", Histogram(probs, 10, 0, 1))
	//normHist := GetNormalizedHisto(probs)
	//t.Logf("normHist:%+v", normHist)

	probs := []float32{0.90483475, 0.2512588, 0.8995428}
	t.Logf("Histogram:%+v", Histogram(probs, 10, 0, 1))
	normHist := GetNormalizedHisto(probs)
	t.Logf("normHist:%+v", normHist)
}

func TestSlidingWindowViewArray(t *testing.T) {
	data := []float64{0.1, 0.2, 0.3, 0.4, 0.5}
	result := SlidingWindowViewArray(data, 3)
	t.Logf("result:%+v", result)
}

func TestSlidingWindowView2DArray(t *testing.T) {
	data := [][]float32{
		{0.1, 0.1, 0.1, 0.1, 0.1},
		{0.2, 0.2, 0.2, 0.2, 0.2},
		{0.3, 0.3, 0.3, 0.3, 0.3},
	}
	result := SlidingWindowView2DArray(data, 2, 0)
	t.Logf("result:%+v", result)
}
