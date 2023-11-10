// hdrpng.cpp
// Copyright (C) 2019-2022 Binomial LLC. All Rights Reserved.
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#if _MSC_VER
// For sprintf(), strcpy() 
#define _CRT_SECURE_NO_WARNINGS (1)
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <vector>
#include <cmath>

#include "helpers.h"

using namespace buminiz;

//--------------------------------------------------------------------------------------------------------------------------

using namespace basisu;

int print_usage()
{
	printf("Usage:\n");
	printf("hdrpng p infile.exr outfile.png - Packs HDR .EXR file to HDR .PNG\n");
	printf("hdrpng u infile.png outfile.exr - Unpacks HDR PNG to HDR .EXR file\n");
	printf("hdrpng c file1.exr file2.exr - Compares two .EXR files for equality (as HALF floats, ignoring any NaN/Inf's)\n");
	printf("Note: The TinyEXR library is used to load .EXR files. Not all compression modes are supported by this library.\n");
	return EXIT_FAILURE;
}

int pack_hdr_png(const char* pIn_filename, const char* pOut_filename)
{
	int num_chans = 0;
	
	imagef hdr_test_img;
	bool status = read_exr(pIn_filename, hdr_test_img, num_chans);
	if (!status)
	{
		fprintf(stderr, "Failed reading source EXR file %s\n", pIn_filename);
		return EXIT_FAILURE;
	}

	const uint32_t width = hdr_test_img.get_width(), height = hdr_test_img.get_height();

	printf(".EXR image file dimensions: %ux%u num_chans: %i\n", width, height, num_chans);

	if (!hdr_test_img.clean_pixels())
	{
		fprintf(stderr, "Warning: one or more source pixels was either out of half float range or was NaN/Inf\n");
	}

	float lowest_nonzero_val = 1e+30f, lowest_val = 1e+30f, highest_val = -1e+30f;

	for (uint32_t y = 0; y < hdr_test_img.get_height(); y++)
	{
		for (uint32_t x = 0; x < hdr_test_img.get_width(); x++)
		{
			const vec4F& c = hdr_test_img(x, y);

			for (uint32_t i = 0; i < 3; i++)
			{
				lowest_val = basisu::minimum(lowest_val, c[i]);

				if (c[i] != 0.0f)
					lowest_nonzero_val = basisu::minimum(lowest_nonzero_val, c[i]);

				highest_val = basisu::maximum(highest_val, c[i]);
			}
		}
	}

	printf("Lowest value: %f, lowest non-zero value: %f, highest value: %f, dynamic range: %f\n", lowest_val, lowest_nonzero_val, highest_val, highest_val / lowest_nonzero_val);
	
	if (lowest_val < 0.0f)
		printf("WARNING: Some pixel(s) were negative, these will pack losslessly but may look odd if their magnitudes are too high\n");
		
	uint16_vec orig_half_img(width * 3 * height);
	uint16_vec half_img(width * 3 * height);
	
	int max_shift = 32;

	for (uint32_t y = 0; y < height; y++)
	{
		for (uint32_t x = 0; x < width; x++)
		{
			const vec4F& p = hdr_test_img(x, y);

			for (uint32_t i = 0; i < 3; i++)
			{
				uint32_t h = float_to_half(p[i]);
				uint32_t orig_h = h;

				orig_half_img[(x + y * width) * 3 + i] = (uint16_t)h;
			
				// Rotate sign bit into LSB
				h = rot_left16((uint16_t)h, 1);
				assert(rot_right16((uint16_t)h, 1) == orig_h);
				
				half_img[(x + y * width) * 3 + i] = (uint16_t)h;

				// Determine # of leading zero bits, ignoring the sign bit
				if (h)
				{
					int lz = clz(h) - 16;
					assert(lz >= 0 && lz <= 16);

					assert((h << lz) <= 0xFFFF);

					max_shift = basisu::minimum<int>(max_shift, lz);
				}
			} // i
		} // x
	} // y

	printf("Max leading zeros: %i\n", max_shift);

	uint32_t high_hist[256];
	clear_obj(high_hist);
	
	for (uint32_t y = 0; y < height; y++)
	{
		for (uint32_t x = 0; x < width; x++)
		{
			for (uint32_t i = 0; i < 3; i++)
			{
				uint16_t &hf = half_img[(x + y * width) * 3 + i];
								
				assert(((uint32_t)hf << max_shift) <= 65535);

				hf <<= max_shift;
																
				uint32_t h = (uint8_t)(hf >> 8);
				high_hist[h]++;
			}
		} // x
	} // y

	uint32_t total_vals_used = 0;
	int remap_old_to_new[256];
	for (uint32_t i = 0; i < 256; i++)
		remap_old_to_new[i] = -1;

	for (uint32_t i = 0; i < 256; i++)
	{
		if (high_hist[i] != 0)
		{
			remap_old_to_new[i] = total_vals_used;
			total_vals_used++;
		}
	}

	assert(total_vals_used >= 1);

	printf("Total used high byte values: %u, unused: %u\n", total_vals_used, 256 - total_vals_used);

	bool val_used[256];
	clear_obj(val_used);

	int remap_new_to_old[256];
	for (uint32_t i = 0; i < 256; i++)
		remap_new_to_old[i] = -1;

	int prev_c = -1;
	for (uint32_t i = 0; i < 256; i++)
	{
		if (remap_old_to_new[i] >= 0)
		{
			int c;
			if (total_vals_used <= 1)
				c = remap_old_to_new[i];
			else
			{
				c = (remap_old_to_new[i] * 255 + ((total_vals_used - 1) / 2)) / (total_vals_used - 1);

				assert(c > prev_c);
			}

			assert(!val_used[c]);

			remap_new_to_old[c] = i;

			remap_old_to_new[i] = c;
			prev_c = c;

			printf("%u ", c);

			val_used[c] = true;
		}
	} // i
	printf("\n");

	for (uint32_t y = 0; y < height; y++)
	{
		for (uint32_t x = 0; x < width; x++)
		{
			for (uint32_t c = 0; c < 3; c++)
			{
				uint16_t& v16 = half_img[(x + y * width) * 3 + c];

				uint32_t hb = v16 >> 8;
				uint32_t lb = v16 & 0xFF;
								
				assert(remap_old_to_new[hb] != -1);
				assert(remap_old_to_new[hb] <= 255);
				assert(remap_new_to_old[remap_old_to_new[hb]] == (int)hb);

				hb = remap_old_to_new[hb];
				
				v16 = (uint16_t)((hb << 8) | lb);

				// Validate the unpack procedure that a half-float aware decoder will be doing
				{
					uint32_t rh = remap_new_to_old[v16 >> 8];
					
					uint32_t recovered_hf = (rh << 8) | (v16 & 0xFF);

					recovered_hf >>= max_shift;
					recovered_hf = rot_right16((uint16_t)recovered_hf, 1);

					if (recovered_hf != orig_half_img[(x + y * width) * 3 + c])
					{
						fprintf(stderr, "Unpack validation failed!\n");
						return EXIT_FAILURE;
					}
				}
			}
		} // x
	} // y

	png_hdr_chunk hdr_chunk;
	hdr_chunk.m_shift_amount = (uint8_t)max_shift;
	for (uint32_t i = 0; i < 256; i++)
		hdr_chunk.m_remap_table[i] = (uint8_t)basisu::maximum(0, remap_new_to_old[i]);

	hdr_chunk.init();

	status = save_rgb16_png(pOut_filename, width, height, half_img, hdr_chunk);
	if (!status)
	{
		fprintf(stderr, "Failed writing output PNG!\n");
		return EXIT_FAILURE;
	}

	printf("Wrote 16-bit HDR PNG file %s\n", pOut_filename);
	printf("Success\n");

	return EXIT_SUCCESS;
}

int unpack_hdr_png(const char* pIn_filename, const char* pOut_filename)
{
	uint32_t width = 0, height = 0;
	uint16_vec img16;
	png_hdr_chunk hdr_chunk;
	clear_obj(hdr_chunk);

	if (!load_rgb16_png(pIn_filename, width, height, img16, hdr_chunk))
	{
		fprintf(stderr, "Failed loading PNG file %s\n", pIn_filename);
		return EXIT_FAILURE;
	}

	if (!hdr_chunk.m_png_len)
	{
		fprintf(stderr, "Could not find hdRa chunk - not a HDR PNG file\n");
		return EXIT_FAILURE;
	}

	imagef recovered_img(width, height);

	for (uint32_t y = 0; y < height; y++)
	{
		for (uint32_t x = 0; x < width; x++)
		{
			for (uint32_t c = 0; c < 3; c++)
			{
				uint32_t v = img16[(x + y * width) * 3 + c];

				uint32_t h = v >> 8;
				uint32_t l = v & 0xFF;
								
				h = hdr_chunk.m_remap_table[h];
				
				v = (h << 8) | l;

				v >>= hdr_chunk.m_shift_amount;
				v = rot_right16((uint16_t)v, 1);

				recovered_img(x, y)[c] = half_to_float((half_float)v);
			}

			recovered_img(x, y)[3] = 1.0f;
		}
	}

	if (!write_exr(pOut_filename, recovered_img, 3, 0))
	{
		fprintf(stderr, "Failed writing output filename %s\n", pOut_filename);
		return EXIT_FAILURE;
	}

	printf("Wrote output EXR file %s\n", pOut_filename);
	printf("Success\n");

	return EXIT_SUCCESS;
}

int compare_exrs(const char* pFilename1, const char* pFilename2)
{
	int num_chans1 = 0;
	imagef hdr_img1;
	bool status = read_exr(pFilename1, hdr_img1, num_chans1);
	if (!status)
	{
		fprintf(stderr, "Failed reading source EXR file %s\n", pFilename1);
		return EXIT_FAILURE;
	}

	printf("Loaded EXR file %s, %ux%u\n", pFilename1, hdr_img1.get_width(), hdr_img1.get_height());

	int num_chans2 = 0;
	imagef hdr_img2;
	status = read_exr(pFilename2, hdr_img2, num_chans2);
	if (!status)
	{
		fprintf(stderr, "Failed reading source EXR file %s\n", pFilename2);
		return EXIT_FAILURE;
	}

	printf("Loaded EXR file %s, %ux%u\n", pFilename2, hdr_img2.get_width(), hdr_img2.get_height());

	if ( (hdr_img1.get_width() != hdr_img2.get_width()) ||
		 (hdr_img1.get_height() != hdr_img2.get_height()) )
	{
		fprintf(stderr, "Dimensions differ!\n");
		return EXIT_FAILURE;
	}

	uint32_t a_invalid_count = 0, b_invalid_count = 0, mismatch_count = 0;
	uint32_t a_clamp_count = 0, b_clamp_count = 0;

	printf("Comparing RGB pixels of each image converted to valid HALF floats:\n");

	for (uint32_t y = 0; y < hdr_img1.get_height(); y++)
	{
		for (uint32_t x = 0; x < hdr_img1.get_width(); x++)
		{
			for (uint32_t c = 0; c < 3; c++)
			{
				float a = hdr_img1(x, y)[c];
				float b = hdr_img2(x, y)[c];
								
				if (std::isnan(a) || std::isinf(a))
				{
					a_invalid_count++;
					continue;
				}

				if (std::isnan(b) || std::isinf(b))
				{
					b_invalid_count++;
					continue;
				}

				// Clamp to MAX_HALF_FLOAT before comparison - there's nothing we can do, HDR PNG is a half float format.
				if (a > MAX_HALF_FLOAT)
				{
					a_clamp_count++;
					a = MAX_HALF_FLOAT;
				}

				if (b > MAX_HALF_FLOAT)
				{
					b_clamp_count++;
					b = MAX_HALF_FLOAT;
				}

				uint32_t ha = float_to_half(a);
				uint32_t hb = float_to_half(b);

				if (ha != hb)
				{
					mismatch_count++;
				}
			}
		} // x
	} // y

	printf("Total first image NaN/Inf: %u\n", a_invalid_count);
	printf("Total first image MAX_HALF_FLOAT clamps: %u\n", a_clamp_count);

	printf("Total second image NaN/Inf: %u\n", b_invalid_count);
	printf("Total second image MAX_HALF_FLOAT clamps: %u\n", b_clamp_count);

	printf("Total mismatches: %u\n", mismatch_count);

	if ((mismatch_count == 0) && !a_invalid_count && !b_invalid_count)
	{
		printf("Image RGB values are equal, and there were no NaN's/Inf's in either image.\n");
	}
	else if ((mismatch_count == 0) && (a_invalid_count || b_invalid_count))
	{
		printf("Image RGB values are equal, however one or more NaN's/Inf's were encountered.\n");
	}
	else
	{
		fprintf(stderr, "Image RGB values are NOT equal, comparison failed!\n");
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

int main(int argc, const char** argv)
{
#ifdef _DEBUG
	printf("DEBUG\n");
#endif

	if (argc != 4)
		return print_usage();

	int mode = tolower(argv[1][0]);
	if ((mode != 'p') && (mode != 'u') && (mode != 'c'))
		return print_usage();

	const char* pIn_filename = argv[2];
	const char* pOut_filename = argv[3];

	if (mode == 'p')
		return pack_hdr_png(pIn_filename, pOut_filename);
	else if (mode == 'u')
		return unpack_hdr_png(pIn_filename, pOut_filename);
	else if (mode == 'c')
		return compare_exrs(pIn_filename, pOut_filename);
}

