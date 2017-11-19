
#include "stdafx.h"
#include "ncorr.h"
#include <direct.h>

using namespace ncorr;
using namespace std;

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

void exportStrainData(string outputFilesDirectory, string csvExport, strain_analysis_output strain_output) {
	// only reachable before Eulerian if.
	string dataPath = outputFilesDirectory + "data/";
	remove(dataPath.c_str());
	_mkdir(dataPath.c_str());
	vector<string> exportVaribles = split(csvExport, ",", false);
	cout << csvExport << endl;
	for (int idx = 0; idx < exportVaribles.size(); idx++) {

		if (exportVaribles[idx].compare("mean exx") == 0) {
			string exxFileName = dataPath + "exx.txt";
			remove(exxFileName.c_str());
			/*ofstream datafile(exxFileName.c_str());

			if (datafile.is_open()) {
			//cout<<"data file open"<<endl;
			for (int idx2 = 0; idx2 < strain_output.strains.size();
			idx2++) {
			Data2D exx_strains = strain_output.strains[idx2].get_exx();
			datafile << getMean(exx_strains) << "\n";

			}
			}

			datafile.close();*/

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
		if (exportVaribles[idx].compare("mean eyy") == 0) {
			string exxFileName = dataPath + "eyy.txt";
			remove(exxFileName.c_str());
			/*ofstream datafile(exxFileName.c_str());
			if (datafile.is_open()) {
			//cout<<"data file open"<<endl;
			for (int idx2 = 0; idx2 < strain_output.strains.size();
			idx2++) {
			Data2D eyy_strains = strain_output.strains[idx2].get_eyy();
			datafile << getMean(eyy_strains) << "\n";

			}
			}

			datafile.close(); */
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
		if (exportVaribles[idx].compare("mean exy") == 0) {
			string exxFileName = dataPath + "exy.txt";
			remove(exxFileName.c_str());
			/*ofstream datafile(exxFileName.c_str());
			if (datafile.is_open()) {
			//cout<<"data file open"<<endl;
			for (int idx2 = 0; idx2 < strain_output.strains.size();
			idx2++) {

			Data2D exy_strains = strain_output.strains[idx2].get_exy();
			datafile << getMean(exy_strains) << "\n";

			}
			}

			datafile.close();*/

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
		if (exportVaribles[idx].compare("mean e1") == 0) {
			string exxFileName = dataPath + "e1.txt";
			remove(exxFileName.c_str());
			/*ofstream datafile(exxFileName.c_str());
			if (datafile.is_open()) {
			//cout<<"data file open"<<endl;
			for (int idx2 = 0; idx2 < strain_output.strains.size();
			idx2++) {

			Data2D exy_strains = strain_output.strains[idx2].get_e1();
			datafile << getMean(exy_strains) << "\n";

			}
			}

			datafile.close();*/

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
	vector<string> exports = split(videoExport, ",", false);
	for (int idx = 0; idx < exports.size(); idx++) {
		if (exports[idx].compare("v") == 0) {
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

		if (exports[idx].compare("u") == 0) {

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
		if (exports[idx].compare("eyy") == 0) {
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

		if (exports[idx].compare("exy") == 0) {
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
		if (exports[idx].compare("exx") == 0) {

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

		if (exports[idx].compare("e1") == 0) {
			double maxCalcStrain = 0;
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {
				Data2D e1_strains = strain_output.strains[i - 1].get_e1();
				maxCalcStrain = max(maxCalcStrain,
						prctile(e1_strains.get_array(), 0.99));
			}
			for (int i = 1; i < strain_input.DIC_input.imgs.size(); i++) {

				string imageName = getFileNameFromPath(
					getImagesArray[i + 1]);
				string saveImagePath = videoPath + imageName + "_e1_"
					+ strainType + ".jpg";
				Data2D e1_strains = strain_output.strains[i - 1].get_e1();

				double maxE1;
				if (strainMax < 0) {
					maxE1 = maxCalcStrain;
				}
				else {
					maxE1 = strainMax;
				}
				save_ncorr_data_over_img(saveImagePath,
					strain_input.DIC_input.imgs[i - 1], e1_strains,
					0.5, 0, maxE1, true, false, false,
					strain_input.DIC_output.units,
					strain_input.DIC_output.units_per_pixel, 50, 1.0,
					11, cv::COLORMAP_JET);
			}
		}
		std::cout << "Made it here 10" << std::endl;
		if (exports[idx].compare("e2") == 0) {
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
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		throw std::invalid_argument(
			"Must have 2 command line input of either 'calculate' or 'load' and the path of the job file");
	}

	// Initialize DIC and strain information ---------------//
	DIC_analysis_input DIC_input;
	DIC_analysis_output DIC_output;
	strain_analysis_input strain_input;
	strain_analysis_output strain_output;
	ifstream infile(argv[2]);
	std::string str;
	std::string file_contents;
	while (std::getline(infile, str)) {
		file_contents += str;
		file_contents.push_back('\n');
	}
	//cout<< argv[1];

	vector<string> firstLevelSplit = split(file_contents, "\u20AC", true);

	//get name
	vector<string> getNameArray = split(firstLevelSplit[1], "\n", true);
	string nameOfFile = split(getNameArray[0], "\t", true)[1];

	//get images
	std::vector<Image2D> imgs;
	vector<string> getImagesArray = split(firstLevelSplit[2], "\n", false);
	//string directoryInUse = getDirPath(getImagesArray[1]);
	string outputFilesDirectory = split(split(firstLevelSplit[5], "\n", false)[23], "\t")[1];

	//return 0;
	for (int i = 1; i < getImagesArray.size(); i++) {

		imgs.push_back(getImagesArray[i]);
		//cv::imshow("image",cv::imread(getImagesArray[i]));
		// cv::waitKey();
	}
	//get roi
	string roiPath;
	roiPath = split(firstLevelSplit[3], "\n", false)[1];

	std::cout << roiPath << std::endl;

	//get DIC parameters
	double scaleFactor;
	INTERP interpType;
	string interpolationType;
	int numThreads;
	SUBREGION subRegion;
	int radius;
	DIC_analysis_config config_DIC_analysis;

	//Get Strain output parameters
	bool isEulerian;
	string csvExport;
	string videoExport;
	string units;
	double units_per_px;
	int fps;
	string openCV_color;
	double end_delay;
	string fourccInput;
	bool colorbar;
	bool axes;
	bool scalebar;
	int numUnits;
	int fontSize;
	int tickMarks;
	double strainMin;
	double strainMax;
	double dispMin;
	double dispMax;
	int strainRadius;
	SUBREGION subRegionStrain;
	string outputType;

	//Get Load Results
	string loadDIC_inputPath;
	string loadDIC_outputPath;
	string loadStrain_inputPath;
	string loadStrain_outputPath;

	scaleFactor = atof(
		split(split(firstLevelSplit[4], "\n", false)[1], "\t")[1].c_str());

	interpolationType =
		split(split(firstLevelSplit[4], "\n", false)[2], "\t")[1];

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

	numThreads = atof(
		split(split(firstLevelSplit[4], "\n", false)[3], "\t")[1].c_str());
	string subRegionstr;
	subRegionstr = split(split(firstLevelSplit[4], "\n", false)[4], "\t")[1];

	if (subRegionstr.compare("Circle") == 0) {
		subRegion = SUBREGION::CIRCLE;
	}
	if (subRegionstr.compare("Nearest") == 0) {
		subRegion = SUBREGION::SQUARE;
	}

	radius = atof(
		split(split(firstLevelSplit[4], "\n", false)[5], "\t")[1].c_str());

	string configurationString;
	configurationString =
		split(split(firstLevelSplit[4], "\n", false)[6], "\t")[1];
	if (configurationString.compare("NO_UPDATE") == 0) {
		config_DIC_analysis = DIC_analysis_config::NO_UPDATE;
	}
	else if (configurationString.compare("KEEP_MOST_POINTS") == 0) {
		config_DIC_analysis = DIC_analysis_config::KEEP_MOST_POINTS;
	}
	else{
		config_DIC_analysis = DIC_analysis_config::REMOVE_BAD_POINTS;
	}

	//Output settings
	isEulerian =
		!(split(split(firstLevelSplit[5], "\n", false)[1], "\t")[1].compare(
			"Lagrangian") == 0);

	csvExport = split(split(firstLevelSplit[5], "\n", false)[2], "\t")[1];

	videoExport = split(split(firstLevelSplit[5], "\n", false)[3], "\t")[1];

	units = split(split(firstLevelSplit[5], "\n", false)[4], "\t")[1];

	units_per_px = atof(
		split(split(firstLevelSplit[5], "\n", false)[5], "\t")[1].c_str());

	fps = atof(
		split(split(firstLevelSplit[5], "\n", false)[6], "\t")[1].c_str());

	openCV_color = split(split(firstLevelSplit[5], "\n", false)[7], "\t")[1];

	end_delay = atof(
		split(split(firstLevelSplit[5], "\n", false)[8], "\t")[1].c_str());

	fourccInput = split(split(firstLevelSplit[5], "\n", false)[9], "\t")[1];

	colorbar =
		(split(split(firstLevelSplit[5], "\n", false)[10], "\t")[1].compare(
			"true") == 0);

	axes = (split(split(firstLevelSplit[5], "\n", false)[11], "\t")[1].compare(
		"true") == 0);

	scalebar =
		(split(split(firstLevelSplit[5], "\n", false)[12], "\t")[1].compare(
			"true") == 0);

	numUnits = atof(
		split(split(firstLevelSplit[5], "\n", false)[13], "\t")[1].c_str());

	fontSize = atof(
		split(split(firstLevelSplit[5], "\n", false)[14], "\t")[1].c_str());

	tickMarks = atof(
		split(split(firstLevelSplit[5], "\n", false)[15], "\t")[1].c_str());

	strainMin = atof(
		split(split(firstLevelSplit[5], "\n", false)[16], "\t")[1].c_str());

	strainMax = atof(
		split(split(firstLevelSplit[5], "\n", false)[17], "\t")[1].c_str());

	dispMin = atof(
		split(split(firstLevelSplit[5], "\n", false)[18], "\t")[1].c_str());

	dispMax = atof(
		split(split(firstLevelSplit[5], "\n", false)[19], "\t")[1].c_str());

	strainRadius = atof(
		split(split(firstLevelSplit[5], "\n", false)[20], "\t")[1].c_str());

	outputType = split(split(firstLevelSplit[5], "\n", false)[22], "\t")[1];

	//return 0;
	string subRegionstrStrain;
	subRegionstrStrain =
		split(split(firstLevelSplit[5], "\n", false)[21], "\t")[1];

	if (subRegionstrStrain.compare("Circle") == 0) {
		subRegionStrain = SUBREGION::CIRCLE;
	}
	if (subRegionstrStrain.compare("Square") == 0) {
		subRegionStrain = SUBREGION::SQUARE;
	}

	loadDIC_inputPath =
		split(split(firstLevelSplit[6], "\n", false)[1], "\t")[1];

	loadDIC_outputPath =
		split(split(firstLevelSplit[6], "\n", false)[2], "\t")[1];

	loadStrain_inputPath = split(split(firstLevelSplit[6], "\n", false)[3],
		"\t")[1];

	loadStrain_outputPath = split(split(firstLevelSplit[6], "\n", false)[4],
		"\t")[1];

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
		processImagesFromBinary(DIC_input, DIC_output, strain_input, strain_output, loadDIC_inputPath, loadDIC_outputPath, strainRadius, subRegionStrain,  interpType,  units,  units_per_px);
	}
	else if (input == "calculate") {

		// Set DIC_input
		std::chrono::time_point<std::chrono::system_clock> start_setUp = std::chrono::system_clock::now();
		std::cout << roiPath << std::endl;
		//cv::imshow("roi", cv::imread(imgs[0].get_filename(), 0));
		//cv::waitKey(0);
		//std::cout << "OpenCV " << cv::getBuildInformation() << endl;
		std::cout << scaleFactor << " " << interpolationType << " " << subRegionstr << " " <<
			numThreads << " " << configurationString << std::endl;
		DIC_input = DIC_analysis_input(imgs, 						// Images
			ROI2D(Image2D(roiPath).get_gs() > 0.5),		// ROI
			scaleFactor,                                      // scalefactor
			interpType,			// Interpolation
			subRegion,					// Subregion shape
			radius,                                      // Subregion radius
			numThreads,                                      // # of threads
			config_DIC_analysis, // DIC configuration for reference image updates
			false);							// Debugging enabled/disabled

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

		if (isEulerian) {
			processImagesFromBinary(DIC_input_eulerian, DIC_output_eulerian, strain_input_eulerian, strain_output_eulerian, saveDIC_inputPath, saveDIC_outputPath, strainRadius, subRegionStrain, interpType, units, units_per_px);
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
	if (isEulerian) {
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
		if(isEulerian)
			saveImages(videoExport, strain_input_eulerian, getImagesArray, videoPath, strainType, DIC_output_eulerian, strain_output_eulerian, strainMax);
		else
			saveImages(videoExport, strain_input, getImagesArray, videoPath, strainType, DIC_output, strain_output, strainMax);

	}
	
	
	return 0;
}