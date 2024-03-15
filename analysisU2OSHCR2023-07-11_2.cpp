#include "ImageAnalysis.hpp"
#include "Nd2ReadSdk.h"


int main() {

	FILE* dots_csv = fopen("../dots_0711_new.csv", "w");
	FILE* nuclei_csv = fopen("../nuclei_0711_new.csv", "w");
	fprintf(dots_csv, "x,y,z,cn_id,cn_x,cn_y,cn_z,channel,condition,pos,intensity\n");
	fprintf(nuclei_csv, "cn_id,cn_x,cn_y,cn_z,zmin,zmax,size,validSize,sum405,sum488,sum561,sum594,condition,pos\n");
	fclose(dots_csv);
	fclose(nuclei_csv);


	std::tuple<const char*, std::vector<int> > condition_stage[9] = {
		std::make_tuple("FISH_shr1x", std::vector<int> {1,2,3,4,5,6,7,8,9,10}),
		std::make_tuple("FISH_shr4x", std::vector<int> {1,2,3,4,5,6,7,8,9,10}),
		std::make_tuple("FISH_0x", std::vector<int> {1,2,3,4,5,6,7,8,9,10}),
	};

	// PARAMETERS
	char* base_dir = "../U2OSFishvsHCR2023-06-29";
	// segmentation parameters
	int threshold_405_higher = 111;
	int threshold_405_lower = 110;
	int sigmaxy = 50;
	int sigmaz = 2;


	// 594 dot parameters
	int threshold_594_lower = 410;
	int threshold_594_higher = 80000;
	int sigmaxy_594 = 1;
	int sigmaz_594 = 1;
	// 647 dot parameters
	int threshold_647_lower_1 = 300;
	int threshold_647_higher_1 = 1200;
	int sigmaxy_647_1 = 30;
	int sigmaz_647_1 = 3;


	int threshold_647_lower_2 = 400;
	int threshold_647_higher_2 = 1200;
	int sigmaxy_647_2 = 5;
	int sigmaz_647_2 = 3;
	// cell dot matching parameters
	int threshold_distSq = 20 * 20;
	int cell_diameter = 200;
	// display parameters 
	float blue_pixel_lower = 100;
	float blue_pixel_upper = 200;
	float green_pixel_lower = 100;
	float green_pixel_upper = 800;
	float red_pixel_lower = 0;
	float red_pixel_upper = 3000;
	float white_pixel_lower = 0;
	float white_pixel_upper = 1;
	int yellow_pixel_lower = 110;
	int yellow_pixel_upper = 3000;
	int gfp_pixel_lower = 130;
	int gfp_pixel_upper = 800;

	for (int i = 0; i < 9; i++) {
		const char* condition;
		std::vector<int> positions;
		std::tie(condition, positions) = condition_stage[i];

		for (int p : positions) {
			char nd2_filename[500];
			sprintf(nd2_filename, "%s\\%s_s%d.nd2", base_dir, condition, p);
			std::cout << "Checking if " << nd2_filename << " exists." << std::endl;


			char mask_filename[500];
			sprintf(mask_filename, "%s\\%s_s%d_640_cp_masks.tif", base_dir, condition, p);

			int width, height, depth;
			width = 2048;
			height = 2048;
			getNd2Depth(nd2_filename, depth);

			printf("Image width: %d, height: %d, depth: %d\n", width, height, depth);

			std::cout << "Loading 561...";
			uint16_t* stack561 = new uint16_t[width * height * depth];
			//complex* stack561fft = new complex[width * height * depth];

			loadNd2(stack561, nd2_filename, width, height, depth, 5, 4);
			std::cout << " done." << std::endl;

			uint16_t* stack594 = new uint16_t[width * height * depth];
			//complex* stack594fft = new complex[width * height * depth];
			//complex* stack594ifft = new complex[width * height * depth];
			float* float594 = new float[width * height * depth];
			uint16_t* median594 = new uint16_t[width * height * depth];
			float* gaussian594 = new float[width * height * depth];
			std::vector<std::tuple<int, int, int>> maxima594;

			std::cout << "Loading 594...";
			// colorindex = 0: Dapi, 1: EGFP, 2: 594, 3: 640, 4:561
			loadNd2(stack594, nd2_filename, width, height, depth, 5, 2);
			std::cout << " done." << std::endl;


			// 647 section
			uint16_t* stack647 = new uint16_t[width * height * depth];
			uint16_t* median647 = new uint16_t[width * height * depth];
			float* gaussian647 = new float[width * height * depth];
			//complex* stack647fft = new complex[width * height * depth];
			//complex* stack647ifft = new complex[width * height * depth];
			float* float647 = new float[width * height * depth];

			std::cout << "Loading 647...";
			loadNd2(stack647, nd2_filename, width, height, depth, 5, 3);
			std::cout << " done." << std::endl;


			std::cout << "Loading masks...";
			uint16_t* masks = new uint16_t[width * height * depth];
			loadTiff(masks, mask_filename, width, height, depth);
			std::cout << "done." << std::endl;

			time_t start, end;

			/*
			std::cout << "performing 561 fft...";
			time(&start);
			fft3D(stack561, {2,2,2,2,2,2,2,2,2,2,2}, { 2,2,2,2,2,2,2,2,2,2,2 }, {2,2,2,3,3}, stack561fft);
			time(&end);
			std::cout << "done." << std::endl;
			std::cout << "FFT took " << end - start << " seconds" << std::endl;
			
			std::cout << "performing 594 fft...";
			time(&start);
			fft3D(stack594, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,3,3 }, stack594fft);
			time(&end);
			std::cout << "done." << std::endl;
			std::cout << "FFT took " << end - start << " seconds" << std::endl;
		
			std::cout << "performing 647 fft...";
			time(&start);
			fft3D(stack647, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,3,3 }, stack647fft);
			time(&end);
			std::cout << "done." << std::endl;
			std::cout << "FFT took " << end - start << " seconds" << std::endl;
			*/
			//std::tuple<int, int, int> shift594 = phase_correlation(stack594fft, stack561fft, width, height, depth);
			//std::tuple<int, int, int> shift594 = { -17,0,0 };
			//shift_fft(stack594fft, shift594, width, height, depth);
			//std::cout << "Checking if phase correlation fixed the shift:" << std::endl;
			//phase_correlation(stack594fft, stack561fft, width, height, depth);

			//std::tuple<int, int, int> shift647 = phase_correlation(stack647fft, stack561fft, width, height, depth);
			//std::tuple<int, int, int> shift647 = { 0,0,0 };
			//shift_fft(stack647fft, shift647, width, height, depth);
			//std::cout << "Checking if phase correlation fixed the shift:" << std::endl;
			//phase_correlation(stack647fft, stack561fft, width, height, depth);

			//std::cout << "performing 594 ifft...";
			//time(&start);
			//fft3D(stack594fft, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,3,3 }, stack594ifft, true);
			//time(&end);
			//std::cout << "done." << std::endl;
			//std::cout << "FFT took " << end - start << " seconds" << std::endl;


			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						int i = width*height*z + width*y + x;
						if (x < width - 17) {
							float594[i] = (float)stack594[i + 17];
							stack594[i] = stack594[i + 17];

						}
						else {
							float594[i] = 0;
							stack594[i] = 0;
						}
					}
				}
			}

			/*
			std::cout << "performing 647 ifft...";
			time(&start);
			fft3D(stack647fft, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,2,2,2,2,2,2,2,2 }, { 2,2,2,3,3 }, stack647ifft, true);
			time(&end);
			std::cout << "done." << std::endl;
			std::cout << "FFT took " << end - start << " seconds" << std::endl;
			*/
			for (int i = 0; i < width * height * depth; i++) {
				//printf("%e + j %e\n", stack594ifft[i].real, stack594ifft[i].imaginary);
				if (stack647[i] > 100) {
					float647[i] = ((float)stack647[i]) - 100;
				}
				else {
					float647[i] = 0;
				}
				if (float594[i] > 100) {
					float594[i] = float594[i] - 100;
				}
				else {
					float594[i] = 0;
				}
				if (stack561[i] > 100) {
					stack561[i] = stack561[i] - 100;
				}
				else {
					stack561[i] = 0;
				}
			}

			std::cout << "Subtracting 561 background from 594 and 647" << std::endl;
			//subtract_projection(float594, stack561, width, height, depth);
			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						int i = width * height * z + width * y + x;
						// float594[i] -= 0.0763308* stack561[i];
						float594[i] -= 0.081 * stack561[i];
					}
				}
			}

			//subtract_projection(float647, stack561, width, height, depth);
			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						int i = width * height * z + width * y + x;
						float647[i] -= 0.0135888 * stack561[i];
					}
				}
			}

			std::cout << "Loading oligoDT...";
			uint16_t* stack = new uint16_t[width * height * depth];
			loadNd2(stack, nd2_filename, width, height, depth, 5, 0);
			std::cout << " done." << std::endl;

			std::vector<Nucleus*> nuclei1;
			uint16_t num_masks = 0;

			float* fpointer = float647;
			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						if ((*fpointer) <= threshold_647_lower_1) {
							*fpointer = 0;
						}
						else if ((*fpointer) >= threshold_647_higher_1) {
							*fpointer = threshold_647_higher_1 - threshold_647_lower_1;
						}
						else {
							*fpointer = (*fpointer) - threshold_647_lower_1;
						}
						fpointer++;
					}
				}
			}

			std::cout << "Converting masks..." << std::endl;
			for (int i = 0; i < width * depth * height; i++) {
				if (masks[i] > num_masks) {
					num_masks = masks[i];
				}
			}
			std::cout << "Found " << num_masks << " masks." << std::endl;

			nuclei1.resize(num_masks);
			for (int i = 0; i < num_masks; i++) {
				Nucleus* blob = new Nucleus(1e12, 10);
				blob->id = i + 1;
				blob->points.reserve(2000000);
				blob->boundary.reserve(200000);
				blob->local_max = std::make_tuple(0, 0, 0);

				nuclei1.at(i) = blob;
			}

			std::cout << "Digesting masks... ";
			time(&start);
			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						int mask = masks[x + width * y + width * height * z];
						if (mask != 0) {
							
							if (x == width - 1 || x == 0 || y == height - 1 || y == 0 || z == depth - 1 || z == 0) {
								nuclei1.at(mask - 1)->boundary.push_back(std::make_tuple(x, y, z));
							}
							else if (
								(masks[x + 1 + width * y + width * height * z] == mask) &&
								(masks[x - 1 + width * y + width * height * z] == mask) &&
								(masks[x + width * (y + 1) + width * height * z] == mask) &&
								(masks[x + width * (y - 1) + width * height * z] == mask) &&
								(masks[x + width * y + width * height * (z + 1)] == mask) &&
								(masks[x + width * y + width * height * (z - 1)] == mask)
								) {
								nuclei1.at(mask - 1)->points.push_back(std::make_tuple(x, y, z));
							}
							else {
								nuclei1.at(mask - 1)->boundary.push_back(std::make_tuple(x, y, z));
							}
						}
					}
				}
			}
			std::cout << "Calculating local max" << std::endl; 
			for (int i = 0; i < num_masks; i++) {
				Nucleus* nuc = nuclei1.at(i);
				float m_x=0, m_y=0, m_z=0;
				for (int j = 0; j < nuc->points.size(); j++) {
					int x, y, z;
					std::tie(x, y, z) = nuc->points.at(j);
					m_x += x;
					m_y += y;
					m_z += z;
				}
				for (int j = 0; j < nuc->boundary.size(); j++) {
					int x, y, z;
					std::tie(x, y, z) = nuc->boundary.at(j);
					m_x += x;
					m_y += y;
					m_z += z;
				}
				float total_points = nuc->points.size() + nuc->boundary.size();
				nuclei1.at(i)->local_max = std::make_tuple((int) (m_x / total_points), (int) (m_y / total_points), (int) (m_z / total_points));
			}

			std::vector<Nucleus*> temp;
			temp.swap(nuclei1);
			for (int i = 0; i < temp.size(); i++) {
				if (temp.at(i)->validSize()) {
					nuclei1.push_back(temp.at(i));
					temp.at(i) = 0;
				}
			}
			time(&end);

			std::cout << "done." << std::endl;
			std::cout << "digesting blob segmentation took: " << difftime(end, start) << " seconds." << std::endl;
			std::cout << nuclei1.size() << " final blobs" << std::endl;
		



			fpointer = float594;
			for (int z = 0; z < depth; z++) {
				for (int y = 0; y < height; y++) {
					for (int x = 0; x < width; x++) {
						if ((*fpointer) <= threshold_594_lower) {
							*fpointer = 0;
						}
						else if ((*fpointer) > threshold_594_higher) {
							*fpointer = threshold_594_higher - threshold_594_lower;
						}
						else {
							*fpointer = (*fpointer) - threshold_594_lower;
						}
						fpointer++;
					}
				}
			}
			std::cout << "computing 594 gaussian...";
			//gaussian_filter3D_parallel(float594, width, height, depth, sigmaxy_594, sigmaz_594, gaussian594);
			std::cout << "done." << std::endl;

			findMaxima(float594, width, height, depth, maxima594);
			std::cout << "Found " << maxima594.size() << " dots." << std::endl;

			if (nuclei1.size() > 0) {
				std::cout << "Matching 594 dots to cells...";
				time(&start);
				for (int i = 0; i < maxima594.size(); i++) {
					int x, y, z;
					std::tie(x, y, z) = maxima594.at(i);
					int mask = masks[x + width * y + width * height * z];
					for (int j = 0; j < nuclei1.size(); j++) {
						if (nuclei1.at(j)->id == mask) {
							nuclei1.at(j)->close_points594.push_back(maxima594.at(i));
						}
					}
				}
				time(&end);
				std::cout << "done. Took " << end - start << " seconds" << std::endl;
			}

			std::cout << "Loading 488...";
			uint16_t* stack488 = new uint16_t[width * height * depth];
			loadNd2(stack488, nd2_filename, width, height, depth, 5, 1); 
			std::cout << " done." << std::endl;

			//std::vector<Nucleus*> nuclei_high;
			//std::vector<int> high_idx;

			std::cout << "Writing data to csv files...";
			FILE* dots_csv = fopen("../dots_0711_new.csv", "a");
			FILE* nuclei_csv = fopen("../nuclei_0711_new.csv", "a");
			for (int i = 0; i < nuclei1.size(); i++) {
				int cn_x, cn_y, cn_z;
				Nucleus nuc = *nuclei1.at(i);
				std::tie(cn_x, cn_y, cn_z) = nuc.local_max;
				double sum488 = 0;
				double sum561 = 0;
				double sum405 = 0;
				double sum594 = 0;
				for (int j = 0; j < nuc.points.size(); j++) {
					int x, y, z;
					std::tie(x, y, z) = nuc.points.at(j);
					sum488 += stack488[x + width * y + width * height * z];
					sum561 += stack561[x + width * y + width * height * z];
					sum405 += stack[x + width * y + width * height * z];
					sum594 += stack594[x + width * y + width * height * z];

				}

				//if (sum405 / nuc.points.size() > 1000) {
					//nuclei_high.push_back(nuclei1.at(i));
					//high_idx.push_back(i);
				//}

				int zmin = cn_z;
				int zmax = cn_z;
				for (int p = 0; p < nuc.boundary.size(); p++) {
					int this_z = std::get<2>(nuc.boundary.at(p));
					if (this_z > zmax) {
						zmax = this_z;
					}
					if (this_z < zmin) {
						zmin = this_z;
					}
				}

				fprintf(nuclei_csv, "%d,%d,%d,%d,%d,%d,%d,%d,%e,%e,%e,%e,%s,%d\n", nuc.id, cn_x, cn_y, cn_z, zmin, zmax, nuc.points.size(), nuc.validSize(), sum405, sum488, sum561, sum594, condition, p);


				for (int j = 0; j < nuc.close_points594.size(); j++) {
					int x, y, z;
					std::tie(x, y, z) = nuc.close_points594.at(j);
					float intensity = float594[width * height * z + width * y + x];
					fprintf(dots_csv, "%d,%d,%d,%d,%d,%d,%d,%d,%s,%d,%e\n", x, y, z, nuc.id, cn_x, cn_y, cn_z, 594, condition, p, intensity);
				}

			}
			fclose(dots_csv);
			fclose(nuclei_csv);
			std::cout << "Done." << std::endl;

			//displayBlobsInt(stack, width, height, depth, nuclei);
			
			

			/*
			  std::cout << "Displaying " << nuclei1.size() << " cells" << std::endl; 			
			int k = 0; 
			int key = 0;
			i = 0;
			do {
				int nx, ny, nz;
				Nucleus* nuc = nuclei1.at(k);
				std::tie(nx, ny, nz) = nuc->local_max;
				cv::Mat img(401, 401, CV_32FC3, cv::Scalar(0, 0, 0));

				// we are going to look at range from maxima-150 to maxima+150
				for (int rx = -200; rx <= 200; rx++) {
					if (nx + rx < 0 || nx + rx >= width) {
						continue;
					}
					for (int ry = -200; ry <= 200; ry++) {
						if (ny + ry < 0 || ny + ry >= height) {
							continue;
						}
						BGR_float & bgr = img.ptr<BGR_float>(ry + 200)[rx + 200];

						float green = float594[(nx + rx) + (ny + ry) * width + width * height * i];
						if (green < green_pixel_lower) bgr.green = 0;
						else if (green > green_pixel_upper) bgr.green = 1;
						else bgr.green = (green - green_pixel_lower) * (1.0 / green_pixel_upper);

						float blue = stack[(nx + rx) + (ny + ry) * width + width * height * i];
						if (blue < blue_pixel_lower) bgr.blue = 0;
						else if (blue > blue_pixel_upper) bgr.blue = 1;
						else bgr.blue = (blue - blue_pixel_lower) * (1.0 / blue_pixel_upper);

						float red = float647[(nx + rx) + (ny + ry) * width + width * height * i];
						if (red < red_pixel_lower) bgr.red = 0;
						else if (red > red_pixel_upper) bgr.red = 1;
						else bgr.red = (red - red_pixel_lower) * (1.0 / red_pixel_upper);
					}
				}


				for (int j = 0; j < nuc->boundary.size(); j++) {
					int x, y, z;

					std::tie(x, y, z) = nuc->boundary.at(j);

					if (x - nx > 200 || x - nx < -200 || y - ny > 200 || y - ny < -200) continue;
					if (z == i) {
						BGR_float & bgr = img.ptr<BGR_float>(y - ny + 200)[x - nx + 200];
						bgr.red = 1;
						bgr.green = 1;
						bgr.blue = 1;

					}
				}

				// for all the dots, show their outline
				for (int j = 0; j < nuc->close_points594.size(); j++) {
					int dx, dy, dz;
					std::tie(dx, dy, dz) = nuc->close_points594.at(j);

					if (((dz - i) < 5 && i - dz < 5)) {
						// if we are close enough
						cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 5, cv::Scalar(0, 1, 0), 1, 8, 0);
						cv::circle(img, cv::Point(dx - nx + 200, dy - ny + 200), 6, cv::Scalar(1, 1, 1), 1, 8, 0);

					}

				}



				cv::imshow("Display window", img);
				key = cv::waitKey(0);
				if (key == 's') {
					i = (i + 1) % depth;
				}
				else if (key == 'w') {
					i = (depth + i - 1) % depth;
				}
				else if (key == 'a') {
					do {
						k = (k + 1) % nuclei1.size();
					} while (!nuclei1.at(k)->validSize());
					i = std::get<2>(nuclei1.at(k)->local_max);
				}
				else if (key == 'd') {
					do {
						k = (nuclei1.size() + k - 1) % nuclei1.size();
					} while (!nuclei1.at(k)->validSize());
					i = std::get<2>(nuclei1.at(k)->local_max);
				}
				else {
					//blink = !blink;
				}
			} while (key != 27); */


			std::cout << "Writing maxprojects...";
			char maxproject_dirname[500];
			sprintf(maxproject_dirname, "%s\\maxprojects_new", base_dir);
			//auto created_new_directory
				//= std::filesystem::create_directory(maxproject_dirname);
			for (int k = 0; k < nuclei1.size();k++) {
				if (!nuclei1.at(k)->validSize()) continue;
				//std::cout << "nucleus: " << k << std::endl;
				int cn_x, cn_y, cn_z, cn_id;
				Nucleus* nuc = nuclei1.at(k);
				cn_id = nuc->id;
				cv::Mat img = closeUpCellBlueMaxProjectImg(nuc, stack,
					blue_pixel_lower,
					blue_pixel_upper + 100,
					green_pixel_lower,
					green_pixel_upper,
					red_pixel_lower,
					red_pixel_upper,
					white_pixel_lower,
					white_pixel_upper,
					width, height, depth,
					true);
				char filename[500];
				sprintf(filename, "%s\\%s_s%d_nuc%d_A_blue.png", maxproject_dirname, condition, p, cn_id);
				cv::imwrite(filename, img);


				img = closeUpCellGreenDotMaxProjectImg(nuc, stack594,
					blue_pixel_lower,
					blue_pixel_upper + 100,
					green_pixel_lower,
					green_pixel_upper,
					red_pixel_lower,
					red_pixel_upper,
					white_pixel_lower,
					white_pixel_upper,
					width, height, depth,
					true);
				sprintf(filename, "%s\\%s_s%d_nuc%d_B_GreenDot.png", maxproject_dirname, condition, p, cn_id);
				cv::imwrite(filename, img);

				img = closeUpCellRedBoundaryMaxProjectImg(nuc, stack647,
					blue_pixel_lower,
					blue_pixel_upper + 100,
					green_pixel_lower,
					green_pixel_upper,
					red_pixel_lower,
					red_pixel_upper,
					white_pixel_lower,
					white_pixel_upper,
					width, height, depth,
					true);
				sprintf(filename, "%s\\%s_s%d_nuc%d_C_RedBoundary.png", maxproject_dirname, condition, p, cn_id);
				cv::imwrite(filename, img);



				img = GreenMaxProject(nuc, stack488,
					gfp_pixel_lower,
					gfp_pixel_upper,
					width, height, depth,
					false);
				sprintf(filename, "%s\\%s_s%d_nuc%d_maxproject_488.png", maxproject_dirname, condition, p, cn_id);
				cv::imwrite(filename, img);



				img = YellowMaxProject(nuc, stack561,
					yellow_pixel_lower,
					yellow_pixel_upper,
					width, height, depth,
					false);
				sprintf(filename, "%s\\%s_s%d_nuc%d_maxproject_561.png", maxproject_dirname, condition, p, cn_id);
				cv::imwrite(filename, img);


			}
			std::cout << "Done." << std::endl;




			for (int i = 0; i < nuclei1.size(); i++) {
				delete nuclei1.at(i);
			}
			for (int i = 0; i < temp.size(); i++) {
				if (temp.at(i)) delete temp.at(i);
			}


			delete[] stack;
			delete[] stack488;
			delete[] stack561;
			delete[] float594;
			delete[] float647;
			delete[] stack594;
			delete[] median594;
			delete[] gaussian594;
			delete[] stack647;
			delete[] median647;
			delete[] gaussian647;
			delete[] masks;
		}
	}


}
