package calc

import (
	"math"
)

// Histogram computes the histogram of the given data with specified bins and range.
func Histogram(data []float32, bins int, min, max float32) []int {
	hist := make([]int, bins)
	binWidth := (max - min) / float32(bins)
	epsilon := float32(1e-9) // 微小的偏移量，避免边界值问题
	for _, value := range data {
		if value >= min && value <= max {
			// 防止浮点数精度问题，将偏移量加入计算
			bin := int((value - min + epsilon) / binWidth)
			if bin == bins {
				bin--
			}
			hist[bin]++
		}
	}
	return hist
}

// L2Norm computes the L2 norm of the given data.
func L2Norm(data []int) float64 {
	var sum float64
	for _, value := range data {
		sum += float64(value * value)
	}
	return math.Sqrt(sum)
}

// GetNormalizedHisto computes the normalized histogram of the given probabilities.
func GetNormalizedHisto(probs []float32) []float32 {
	cov := len(probs)
	if cov == 0 {
		return make([]float32, 10)
	}
	hist := Histogram(probs, 10, 0, 1)
	norm := L2Norm(hist)
	normHist := make([]float32, 10)
	for i, value := range hist {
		normHist[i] = float32(value) / float32(norm)
	}
	return normHist
}

// GetNormalizedHistoHiFi computes the normalized histogram of the given probabilities.
func GetNormalizedHistoHiFi(probs []float32) []float32 {
	cov := len(probs)
	if cov == 0 {
		return make([]float32, 30)
	}
	hist := Histogram(probs, 30, 0, 1)
	norm := L2Norm(hist)
	normHist := make([]float32, 30)
	for i, value := range hist {
		normHist[i] = float32(value) / float32(norm)
	}
	return normHist
}

/*// SlidingWindowView creates a sliding window view of the given data with specified window size and axis.
func SlidingWindowView(data *mat.Dense, windowSize int, axis int) []*mat.Dense {
	rows, cols := data.Dims()
	if axis == 0 {
		result := make([]*mat.Dense, rows-windowSize+1)
		for i := 0; i < rows-windowSize+1; i++ {
			result[i] = data.Slice(i, i+windowSize, 0, cols).(*mat.Dense)
		}
		return result
	}
	// For axis == 1
	result := make([]*mat.Dense, cols-windowSize+1)
	for i := 0; i < cols-windowSize+1; i++ {
		result[i] = data.Slice(0, rows, i, i+windowSize).(*mat.Dense)
	}
	return result
}*/

// SlidingWindowView2DArray creates a sliding window view of the given 2D array with specified window size along the given axis.
func SlidingWindowView2DArray(data [][]float32, windowSize int, axis int) [][][]float32 {
	var result [][][]float32
	rows := len(data)
	cols := len(data[0])

	if axis == 0 {
		// Slide along rows
		for i := 0; i <= rows-windowSize; i++ {
			window := make([][]float32, windowSize)
			for j := 0; j < windowSize; j++ {
				window[j] = make([]float32, cols)
				copy(window[j], data[i+j])
			}
			result = append(result, window)
		}
	} else {
		// Slide along columns
		for i := 0; i <= cols-windowSize; i++ {
			window := make([][]float32, rows)
			for j := 0; j < rows; j++ {
				window[j] = make([]float32, windowSize)
				copy(window[j], data[j][i:i+windowSize])
			}
			result = append(result, window)
		}
	}

	return result
}

// SlidingWindowViewIntArray creates a sliding window view of the given int array with specified window size.
func SlidingWindowViewIntArray(data []int, windowSize int) [][]int {
	var result [][]int
	nums := len(data)
	for i := 0; i <= nums-windowSize; i++ {
		window := make([]int, windowSize)
		copy(window, data[i:i+windowSize])
		result = append(result, window)
	}
	return result
}

// SlidingWindowViewArray creates a sliding window view of the given 1D array with the specified window size.
// The function is generic and works with any type T.
func SlidingWindowViewArray[T any](data []T, windowSize int) [][]T {
	var result [][]T
	n := len(data)

	if windowSize > n {
		return result // If the window size is larger than the array, return an empty result
	}

	for i := 0; i <= n-windowSize; i++ {
		window := make([]T, windowSize)
		copy(window, data[i:i+windowSize])
		result = append(result, window)
	}

	return result
}

// Mean Helper function to calculate mean
func Mean(data []float64) float64 {
	sum := 0.0
	for _, val := range data {
		sum += val
	}
	return sum / float64(len(data))
}
