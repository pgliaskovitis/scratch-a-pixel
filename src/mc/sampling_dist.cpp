/*
 * Copyright (C) 2012  www.scratchapixel.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <random>

static const int MAX_NUM = 20; // items in the population are numbered (number between 0 and 20)
static const int MAX_FREQ = 50; // number of items holding a particular number varies between 1 and 50

int main(int argc, char **argv)
{
	if (argc < 3) {
		std::cerr << "Usage: mc_samplingdist <min_samples> <max_samples>" << std::endl;
		return 1;
	}

	int minSampples = atoi(argv[1]); // minimum sample size for each "first level" mean
	int maxSampples = atoi(argv[2]); // maximum sample size for each "first level" mean
	int population[MAX_NUM + 1];
	int meansCount[MAX_NUM + 1];
	int popSize = 0;
	float popMean = 0;
	float popVar = 0;
	float expectedValueDistrMeans = .0f;
	float varianceDistrMeans = .0f;
	static const int numSamples = 10000;
	std::mt19937 rng;

	for (int i = 0; i <= MAX_NUM; ++i) {
		meansCount[i] = 0;
	}

	rng.seed(17);

	// creation population
	std::uniform_int_distribution<uint32_t> distr(1, MAX_FREQ);

	for (int i = 0; i <= MAX_NUM; ++i) {
		population[i] = distr(rng);
		popSize += population[i];
		popMean += population[i] * i; // prob * x_i
		popVar += population[i] * i * i; // prob * x_i^2
	}

	popMean /= popSize;
	popVar /= popSize;
	popVar -= popMean * popMean;
	fprintf(stderr, "size %d mean %f var %f\n", popSize, popMean, popVar);
	std::uniform_int_distribution<uint32_t> n_samples_distr(minSampples, maxSampples);
	std::uniform_int_distribution<uint32_t> pick_item_distr(0, popSize - 1);

	// now that we have some data and stats to work with sample it
	for (int i = 0; i < numSamples; ++i) {
		// int n = n_samples_distr(rng); // variable sample size ok for sampling distribution
		int n = maxSampples; // fix sample size needed for "second level" mean estimation
		float sample_mean = 0, sample_variance = 0;
		// draw samples from population and compute stats
		for (int j = 0; j < n; ++j) {
			int item_index = pick_item_distr(rng), k;
			for (k = 0; k <= MAX_NUM; ++k) {
				item_index -= population[k];
				if (item_index < 0) break;
			}
			// k is the value we picked up from population,
			// this is the outcome a number between [0:20]
			sample_mean += k;
			sample_variance += k * k;
		}
		// "first level" mean estimation
		sample_mean /= n;
		sample_variance /= n;
		sample_variance -= sample_mean * sample_mean;

		float c1[3] = { 1, 0, 0 };
		float c2[3] = { 0, 0, 1 };
		float t = (n - minSampples) / (float)(maxSampples - minSampples);
		float r = c1[0] * (1 - t) + c2[0] * t;
		float g = c1[1] * (1 - t) + c2[1] * t;
		float b = c1[2] * (1 - t) + c2[2] * t;
		fprintf(stderr, "sample mean: %f sample variance: %f col: %f %f %f;\n", sample_mean, sample_variance, r, g, b);

		// sampling distribution pdf
		meansCount[(int)sample_mean]++;
		for (int i = 0; i <= MAX_NUM; ++i) {
			fprintf(stderr, "%d %d\n", i, meansCount[i]);
		}

		// "second level" mean estimation
		expectedValueDistrMeans += sample_mean;
		varianceDistrMeans += sample_mean * sample_mean;
	}

	expectedValueDistrMeans /= numSamples;
	varianceDistrMeans /= numSamples;
	varianceDistrMeans -= expectedValueDistrMeans * expectedValueDistrMeans;
	fprintf(stderr, "Expected Value of the Mean %f Standard Deviation %f\n", expectedValueDistrMeans, sqrt(varianceDistrMeans));

	return 0;
}