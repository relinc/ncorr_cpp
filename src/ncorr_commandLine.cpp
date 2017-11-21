
#include "stdafx.h"
#include "ncorr.h"
#include <direct.h>
#include "json.hpp"

using namespace ncorr;
using namespace std;
using json = nlohmann::json;


vector<string> split(const string& s, const string& delim,
	const bool keep_empty = true) {
	vector<string> result;
	if (delim.empty()) {
		result.push_back(s);
		return result;
	}
	string::const_iterator substart = s.begin(), subend;
	while (true) {
		subend = search(substart, s.end(), delim.begin(), delim.end());
		string temp(substart, subend);
		if (keep_empty || !temp.empty()) {
			result.push_back(temp);
		}
		if (subend == s.end()) {
			break;
		}
		substart = subend + delim.size();
	}
	return result;
}
string getDirPath(const string & filepath) {
	return filepath.substr(0, filepath.find_last_of("\\/") + 1);

}
string getFileNameFromPath(const string & filepath) {
	string filepathWithExtension = filepath.substr(filepath.find_last_of("\\/") + 1, filepath.length());
	return filepathWithExtension.substr(0, filepathWithExtension.find_last_of("."));

}
double getMean(const Data2D & dataOfInterest) {
	double sumValue = 0;
	int size = 0;
	for (int width = 0; width<dataOfInterest.data_width(); width++) {
		for (int height = 0; height<dataOfInterest.data_height(); height++) {
			if (dataOfInterest.get_roi()(height, width)) {
				sumValue += dataOfInterest.get_array()(height, width);
				size++;
			}
		}

	}
	return sumValue / size;
}

void exportStrainData(string outputFilesDirectory, vector<string> csvExport, strain_analysis_output strain_output) {
	// only reachable before Eulerian if.
	string dataPath = outputFilesDirectory + "data/";
	remove(dataPath.c_str());
	_mkdir(dataPath.c_str());

	for (int idx = 0; idx < csvExport.size(); idx++) {

		if (csvExport[idx].compare("mean exx") == 0) {
			string exxFileName = dataPath + "exx.txt";
			remove(exxFileName.c_str());

			FILE* f = fopen(exxFileName.c_str(), "w");
			if (f == NULL)
			{
				printf("Error opening file!\n");
				exit(1);
			}
			else {
				for (int idx2 = 0; idx2 < strain_output.strains.size();
				idx2++) {
					Data2D exx_strains = strain_output.strains[idx2].get_exx();
					double trueStrain = std::log(1 + getMean(exx_strains));
					fprintf(f, "%f\n", trueStrain);

				}
			}

		}
		if (csvExport[idx].compare("mean eyy") == 0) {
			string exxFileName = dataPath + "eyy.txt";
			remove(exxFileName.c_str());

			FILE* f = fopen(exxFileName.c_str(), "w");
			if (f == NULL)
			{
				printf("Error opening file!\n");
				exit(1);
			}
			else {
				for (int idx2 = 0; idx2 < strain_output.strains.size();
				idx2++) {
					Data2D eyy_strains = strain_output.strains[idx2].get_eyy();
					double trueStrain = std::log(1 + getMean(eyy_strains));
					fprintf(f, "%f\n", trueStrain);

				}
			}

		}
		if (csvExport[idx].compare("mean exy") == 0) {
			string exxFileName = dataPath + "exy.txt";
			remove(exxFileName.c_str());


			FILE* f = fopen(exxFileName.c_str(), "w");
			if (f == NULL)
			{
				printf("Error opening file!\n");
				exit(1);
			}
			else {
				for (int idx2 = 0; idx2 < strain_output.strains.size();
				idx2++) {

					Data2D exy_strains = strain_output.strains[idx2].get_exy();
					double trueStrain = std::log(1 + getMean(exy_strains));
					fprintf(f, "%f\n", trueStrain);

				}
			}

		}
		if (csvExport[idx].compare("mean e1") == 0) {
			string exxFileName = dataPath + "e1.txt";
			remove(exxFileName.c_str());

			FILE* f = fopen(exxFileName.c_str(), "w");
			if (f == NULL)
			{
				printf("Error opening file!\n");
				exit(1);
			}
			else {
				for (int idx2 = 0; idx2 < strain_output.strains.size();
				idx2++) {

					Data2D e1_strains = strain_output.strains[idx2].get_e1();
					double trueStrain = std::log(1 + getMean(e1_strains));
					fprintf(f, "%f\n", trueStrain);

				}
			}

		}
	}
}

void processImagesFromBinary(DIC_analysis_input & DIC_input, DIC_analysis_output & DIC_output, strain_analysis_input& strain_input, 
							 strain_analysis_output& strain_output, string loadDIC_inputPath, 
							 string loadDIC_outputPath, int strainRadius, SUBREGION subRegionStrain, INTERP interpType, string units, double units_per_px) {

	

	// Set strain input


	DIC_input = DIC_analysis_input::load(loadDIC_inputPath);
	DIC_output = DIC_analysis_output::load(loadDIC_outputPath);
	DIC_output = change_perspective(DIC_output, interpType);
	// Set units of DIC_output (provide units/pixel)
	DIC_output = set_units(DIC_output, units, units_per_px);
	// Set strain input
	strain_input = strain_analysis_input(DIC_input, DIC_output,
		subRegionStrain,					// Strain subregion shape
		strainRadius);						// Strain subregion radius

											// Perform strain_analysis

	strain_output = strain_analysis(strain_input);
}

void saveImages(string videoExport, strain_analysis_input &strain_input, vector<string> getImagesArray, string videoPath, string strainType, DIC_analysis_output &DIC_output, strain_analysis_output &strain_output, double strainMax)
{
	// vector<string> exports = split(videoExport, ",", false);
	// for (int idx = 0; idx < exports.size(); idx++) {
		if (videoExport.compare("v") == 0) {
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_v_"
					+ strainType + ".jpg";
				Data2D v_dips = DIC_output.disps[i - 1].get_v();
				double minDisp = min(prctile(v_dips.get_array(), 0.01),
					prctile(v_dips.get_array(), 0.01));
				double maxDisp = max(prctile(v_dips.get_array(), 0.99),
					prctile(v_dips.get_array(), 0.99));
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], v_dips, 0.5,
					minDisp, maxDisp, true, true, true,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);

			}

		}

		if (videoExport.compare("u") == 0) {

			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_u_"
					+ strainType + ".jpg";
				Data2D u_dips = DIC_output.disps[i - 1].get_u();
				double minDisp = min(prctile(u_dips.get_array(), 0.01),
					prctile(u_dips.get_array(), 0.01));
				double maxDisp = max(prctile(u_dips.get_array(), 0.99),
					prctile(u_dips.get_array(), 0.99));

				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], u_dips, 0.5,
					minDisp, maxDisp, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);

			}
		}
		std::cout << "Made it here 6" << std::endl;
		if (videoExport.compare("eyy") == 0) {
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_eyy_"
					+ strainType + ".jpg";
				Data2D eyy_strains = strain_output.strains[i - 1].get_eyy();
				double maxEyy;
				if (strainMax < 0) {

					maxEyy = max(prctile(eyy_strains.get_array(), 0.99),
						prctile(eyy_strains.get_array(), 0.99));
				}
				else {
					maxEyy = strainMax;
				}
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], eyy_strains,
					0.5, 0, maxEyy, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);

			}
		}

		if (videoExport.compare("exy") == 0) {
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_exy_"
					+ strainType + ".jpg";
				Data2D exy_strains = strain_output.strains[i - 1].get_exy();
				double maxExy;
				if (strainMax < 0) {

					maxExy = max(prctile(exy_strains.get_array(), 0.99),
						prctile(exy_strains.get_array(), 0.99));
				}
				else {
					maxExy = strainMax;
				}
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], exy_strains,
					0.5, 0, maxExy, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);

			}
		}
		if (videoExport.compare("exx") == 0) {

			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_exx_"
					+ strainType + ".jpg";
				Data2D exx_strains = strain_output.strains[i - 1].get_exx();
				double maxExx;
				if (strainMax < 0) {

					maxExx = max(prctile(exx_strains.get_array(), 0.99),
						prctile(exx_strains.get_array(), 0.99));
				}
				else {
					maxExx = strainMax;
				}
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], exx_strains,
					0.5, 0, maxExx, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);

			}
		}

		if (videoExport.compare("e1") == 0) {
			double maxCalcStrain = 0;
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {
				Data2D e1_strains = strain_output.strains[i - 1].get_e1();
				maxCalcStrain = max(maxCalcStrain,
						prctile(e1_strains.get_array(), 0.99));
			}
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {
				std::cout<<"here"<<getImagesArray[i + 1]<<std::endl;
				string imageName = getFileNameFromPath(getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_e1_"+ strainType + ".jpg";
				Data2D e1_strains = strain_output.strains[i - 1].get_e1();
				std::cout<<"here1"<<std::endl;
				double maxE1;
				if (strainMax < 0) {
					maxE1 = maxCalcStrain;
				}
				else {
					maxE1 = strainMax;
				}
				std::cout<<"here2"<<std::endl;
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], e1_strains,
					0.5, 0, maxE1, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);
			}
		}
		std::cout << "Made it here 10" << std::endl;
		if (videoExport.compare("e2") == 0) {
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_e2_"
					+ strainType + ".jpg";
				Data2D e2_strains = strain_output.strains[i - 1].get_e2();
				double maxE2;
				if (strainMax < 0) {

					maxE2 = max(prctile(e2_strains.get_array(), 0.99),
						prctile(e2_strains.get_array(), 0.99));
				}
				else {
					maxE2 = strainMax;
				}
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], e2_strains,
					0.5, 0, maxE2, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);

			}
		}

	
}

void parseDICInput(DIC_analysis_input& DIC_input,json inputFile,vector<string>& image_names)
{
	//get DIC parameters
	double scaleFactor;
	INTERP interpType;
	string interpolationType;
	int numThreads;
	SUBREGION subRegion;
	int radius;
	DIC_analysis_config config_DIC_analysis;

	scaleFactor = inputFile.find("dic_settings")->find("scale_factor")->get<double>();

	interpolationType = inputFile.find("dic_settings")->find("interpolation")->get<string>();
	if (interpolationType.compare("Quintic B-spline Precompute") == 0) {
		cout << interpolationType << endl << endl;
		interpType = INTERP::QUINTIC_BSPLINE_PRECOMPUTE;
	}
	if (interpolationType.compare("Quintic B-spline") == 0) {
		interpType = INTERP::QUINTIC_BSPLINE;
	}
	if (interpolationType.compare("Cubic Keys Precompute") == 0) {
		interpType = INTERP::CUBIC_KEYS_PRECOMPUTE;
	}
	if (interpolationType.compare("Cubic Keys") == 0) {
		interpType = INTERP::CUBIC_KEYS;
	}
	if (interpolationType.compare("Linear") == 0) {
		interpType = INTERP::LINEAR;
	}
	if (interpolationType.compare("Nearest") == 0) {
		interpType = INTERP::NEAREST;
	}

	numThreads = inputFile.find("dic_settings")->find("threads")->get<int>();
	string subRegionstr = inputFile.find("dic_settings")->find("subregion")->get<string>();

	if (subRegionstr.compare("Circle") == 0) {
		subRegion = SUBREGION::CIRCLE;
	}
	if (subRegionstr.compare("Nearest") == 0) {
		subRegion = SUBREGION::SQUARE;
	}

	radius = inputFile.find("dic_settings")->find("radius_sub")->get<int>();

	string configurationString=inputFile.find("dic_settings")->find("dic_config")->get<string>();
	if (configurationString.compare("NO_UPDATE") == 0) {
		config_DIC_analysis = DIC_analysis_config::NO_UPDATE;
	}
	else if (configurationString.compare("KEEP_MOST_POINTS") == 0) {
		config_DIC_analysis = DIC_analysis_config::KEEP_MOST_POINTS;
	}
	else{
		config_DIC_analysis = DIC_analysis_config::REMOVE_BAD_POINTS;
	}
	std::vector<Image2D> imgs;
	for (json::iterator it = inputFile.find("images")->begin(); it != inputFile.find("images")->end(); ++it) {
  		image_names.push_back(it->get<string>());
	}
	//string directoryInUse = getDirPath(getImagesArray[1]);
	//return 0;
	for (int i = 1; i < image_names.size(); i++) {

		imgs.push_back(image_names[i]);
	}
	string roiPath = inputFile.find("roi")->get<string>();

	DIC_input = DIC_analysis_input(imgs, 						// Images
			ROI2D(Image2D(roiPath).get_gs() > 0.5),		// ROI
			scaleFactor,                                      // scalefactor
			interpType,			// Interpolation
			subRegion,					// Subregion shape
			radius,                                      // Subregion radius
			numThreads,                                      // # of threads
			config_DIC_analysis, // DIC configuration for reference image updates
			false);	

}

int main(int argc, char *argv[]) {
	
	if (argc != 3) {
		throw std::invalid_argument(
			"Must have 2 command line input of either 'calculate' or 'load' and the path of the job file");
	}
	try{
	// Initialize DIC and strain information ---------------//
	DIC_analysis_input DIC_input;
	string nameOfFile;
	
	std::ifstream i(argv[2]);
	json inputFile;
	i >> inputFile;
	i.close();
	std::vector<std::string> image_names;
	parseDICInput(DIC_input,inputFile,image_names);
	

	nameOfFile=inputFile.find("sample_name")->get<string>();

	DIC_analysis_output DIC_output;

	strain_analysis_input strain_input;
	strain_analysis_output strain_output;
	

	string loadDIC_inputPath=inputFile.find("results")->find("dic_input")->get<string>();
	string loadDIC_outputPath=inputFile.find("results")->find("dic_output")->get<string>();
	string loadStrain_inputPath=inputFile.find("results")->find("strain_input")->get<string>();
	string loadStrain_outputPath=inputFile.find("results")->find("strain_output")->get<string>();
	//Get Load Results
	bool eulerian_out=inputFile.find("output_settings")->find("output_strain_mode")->get<string>().compare("Eulerian")==0;
	bool eulerian_video=inputFile.find("output_settings")->find("video_strain_mode")->get<string>().compare("Eulerian")==0;
	vector<string> csvExport;
	json csvExportJson=*inputFile.find("output_settings")->find("csv_out");
	for (json::iterator it = csvExportJson.begin(); it != csvExportJson.end(); ++it) {
  		csvExport.push_back(it->get<string>());
	}
	string videoExport=inputFile.find("output_settings")->find("video_img_out")->get<string>();
	string units=inputFile.find("output_settings")->find("units")->get<string>();
	double units_per_px=inputFile.find("output_settings")->find("units_per_px")->get<double>();
	int fps = inputFile.find("output_settings")->find("fps")->get<int>();
	string openCV_color =inputFile.find("output_settings")->find("opencv_color")->get<string>();
	double end_delay=inputFile.find("output_settings")->find("end_delay")->get<double>();
	json fourccInputJson = *inputFile.find("output_settings")->find("fourcc");
	vector<string> fourccInput;
	
	for (json::iterator it = fourccInputJson.begin(); it != fourccInputJson.end(); ++it) {
  		csvExport.push_back(it->get<string>());
	} 
	
	bool colorbar=inputFile.find("output_settings")->find("colorbar")->get<bool>();
	bool axes = inputFile.find("output_settings")->find("axes")->get<bool>();
	bool scalebar = inputFile.find("output_settings")->find("scalebar")->get<bool>();
	int numUnits = inputFile.find("output_settings")->find("num_units")->get<int>();
	int fontSize = inputFile.find("output_settings")->find("font_size")->get<int>();
	int tickMarks = inputFile.find("output_settings")->find("tick_marks")->get<int>();
	double strainMin = inputFile.find("output_settings")->find("strain_min")->get<double>();
	double strainMax = inputFile.find("output_settings")->find("strain_max")->get<double>();
	double dispMin = inputFile.find("output_settings")->find("disp_min")->get<double>();
	double dispMax = inputFile.find("output_settings")->find("disp_max")->get<double>();
	int strainRadius = inputFile.find("output_settings")->find("strain_radius")->get<int>();
	SUBREGION subRegionStrain;
	string outputType = inputFile.find("output_settings")->find("output")->get<string>();
	string outputFilesDirectory = inputFile.find("output_settings")->find("output_dir")->get<string>();
	//return 0;
	string subRegionstrStrain = inputFile.find("output_settings")->find("subregion")->get<string>();
	if (subRegionstrStrain.compare("Circle") == 0) {
		subRegionStrain = SUBREGION::CIRCLE;
	}
	if (subRegionstrStrain.compare("Square") == 0) {
		subRegionStrain = SUBREGION::SQUARE;
	}

	strain_analysis_input strain_input_lagrange;
	strain_analysis_output strain_output_lagrange;
	strain_analysis_input strain_input_eulerian;
	strain_analysis_output strain_output_eulerian;
	DIC_analysis_output DIC_output_eulerian;
	DIC_analysis_input DIC_input_eulerian;

	// Determine whether or not to perform calculations or 
	// load data (only load data if analysis has already 
	// been done and saved or else throw an exception).
	std::string input(argv[1]);
	if (input == "load") {
		processImagesFromBinary(DIC_input, 
								DIC_output, 
								strain_input, 
								strain_output, 
								loadDIC_inputPath, 
								loadDIC_outputPath, 
								strainRadius, 
								subRegionStrain,  
								DIC_input.interp_type,  
								units,  
								units_per_px);
	}
	else if (input == "calculate") {

		// Set DIC_input
		std::chrono::time_point<std::chrono::system_clock> start_setUp = std::chrono::system_clock::now();

		std::chrono::time_point<std::chrono::system_clock> end_setUp = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds_setUp = end_setUp - start_setUp;
		std::cout << "Time to setup: " << elapsed_seconds_setUp.count() << "." << std::endl;
		// Perform DIC_analysis    
		DIC_output = DIC_analysis(DIC_input);
		//std::cout << "This happened now." << std::endl;
		// Save outputs as binary
		string savePath = outputFilesDirectory + "save/";
		remove(savePath.c_str());
		_mkdir(savePath.c_str());

		string saveDIC_inputPath = savePath + "DIC_input.bin";
		string saveDIC_outputPath = savePath + "DIC_output.bin";

		save(DIC_input, saveDIC_inputPath);
		save(DIC_output, saveDIC_outputPath);

		DIC_analysis_output DIC_output_original = DIC_output;
		//DIC_output = set_units(DIC_output, units, units_per_px);
		// Set strain input
		strain_input_lagrange = strain_analysis_input(DIC_input, DIC_output_original,
			subRegionStrain,					// Strain subregion shape
			strainRadius);						// Strain subregion radius

												// Perform strain_analysis
		strain_output_lagrange = strain_analysis(strain_input_lagrange);

		//export strain data here, and here only.

		exportStrainData(outputFilesDirectory, csvExport, strain_output_lagrange);

		strain_output = strain_output_lagrange;
		DIC_output=DIC_output_original;
		strain_input = strain_input_lagrange;
		
		DIC_output = set_units(DIC_output, units, units_per_px);
		//export images

		if (!eulerian_out&&eulerian_video) {
			processImagesFromBinary(DIC_input_eulerian, 
									DIC_output_eulerian, 
									strain_input_eulerian, 
									strain_output_eulerian, 
									saveDIC_inputPath, 
									saveDIC_outputPath, 
									strainRadius, 
									subRegionStrain, 
									DIC_input.interp_type, 
									units, 
									units_per_px);
		}
	}
	else {
		throw std::invalid_argument(
			"Input of " + input
			+ " is not recognized. Must be either 'calculate' or 'load'");
	}

	// Create Videos ---------------------------------------//
	// Note that more inputs can be used to modify plots. 
	// If video is not saving correctly, try changing the 
	// input codec using cv::VideoWriter::fourcc(...)). Check 
	// the opencv documentation on video codecs. By default, 
	// ncorr uses cv::VideoWriter::fourcc('M','J','P','G')).
	// Save outputs as binary
	std::cout << "Before Video Start" << std::endl;
	string videoPath = outputFilesDirectory + "video/";
	remove(videoPath.c_str());
	_mkdir(videoPath.c_str());
	
	string strainType = "Lagrangian";
	if (eulerian_video) {
		strainType = "Eulerian";
	}
	
	string saveDICVideo_V_inputPath = videoPath + nameOfFile + "_V_"
		+ strainType + ".avi";
	
	string saveDICVideo_U_outputPath = videoPath + nameOfFile + "_U_"
		+ strainType + ".avi";
	string saveStrainVideo_eyy_inputPath = videoPath + nameOfFile + "_eyy_"
		+ strainType + ".avi";
	string saveStrainVideo_exx_outputPath = videoPath + nameOfFile + "_exx_"
		+ strainType + ".avi";
	
	string saveStrainVideo_exy_outputPath = videoPath + nameOfFile + "_exy_"
		+ strainType + ".avi";
	string saveStrainVideo_e1_outputPath = videoPath + nameOfFile + "_e1_"
		+ strainType + ".avi";
	string saveStrainVideo_e2_outputPath = videoPath + nameOfFile + "_e2_"
		+ strainType + ".avi";
	std::cout << "Made it here 1" << std::endl;
	if (outputType.compare("Video") == 0) {
		vector<string> exports = split(videoExport, ",", false);
		
		for (int idx = 0; idx < exports.size(); idx++) {
			if (exports[idx].compare("v") == 0) {
				save_DIC_video(saveDICVideo_V_inputPath, DIC_input, DIC_output,
					DISP::V, 0.5,		// Alpha
					15);		// FPS
			}
			if (exports[idx].compare("u") == 0) {
				save_DIC_video(saveDICVideo_U_outputPath, DIC_input, DIC_output,
					DISP::U, 0.5,		// Alpha
					15);		// FPS
			}
			if (exports[idx].compare("eyy") == 0) {
				save_strain_video(saveStrainVideo_eyy_inputPath, strain_input,
					strain_output, STRAIN::EYY, 0.5,		// Alpha
					15);		// FPS
			}
			if (exports[idx].compare("exy") == 0) {
				save_strain_video(saveStrainVideo_exy_outputPath, strain_input,
					strain_output, STRAIN::EXY, 0.5,		// Alpha
					15,strainMin,strainMax);		// FPS
			}
			if (exports[idx].compare("exx") == 0) {
				save_strain_video(saveStrainVideo_exx_outputPath, strain_input,
					strain_output, STRAIN::EXX, 0.5,		// Alpha
					15); 		// FPS
			}
			if (exports[idx].compare("e1") == 0) {
				save_strain_video(saveStrainVideo_e1_outputPath, strain_input,
					strain_output, STRAIN::E1, 0.5,		// Alpha
					15); 		// FPS
			}
			if (exports[idx].compare("e2") == 0) {
				save_strain_video(saveStrainVideo_e2_outputPath, strain_input,
					strain_output, STRAIN::E1, 0.5,		// Alpha
					15); 		// FPS
			}
		}
	}

	else {
		if(eulerian_video)
			saveImages(videoExport, 
						strain_input_eulerian, 
						image_names, 
						videoPath, 
						strainType, 
						DIC_output_eulerian, 
						strain_output_eulerian, 
						strainMax);
		else
			saveImages(videoExport, 
						strain_input, 
						image_names, 
						videoPath, 
						strainType, 
						DIC_output, 
						strain_output, 
						strainMax);

	}
	}
	catch(const std::exception& e)
	{
		std::cout<<e.what()<<std::endl;
	}
	
	std::cout<<"exiting ncorr"<<std::endl;
	return 0;
}
