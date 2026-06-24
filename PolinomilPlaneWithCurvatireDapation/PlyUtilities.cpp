
#include "PlyUtilities.h"
#include <vector>
#include <fstream>
#include <iomanip>
#include <thread>
#include <mutex>
#include <chrono>
#include <algorithm>

using namespace std;
PlyUtilities::PlyUtilities(string path, string sampleName, int devisor) {
    this->path = path + "\\" + sampleName;
    this->sampleName = sampleName;
	this->createFolderIfNotexist(); 
    this->g_ProfileIndex = 0;
    this->scanPerMM = devisor / 1000.0;
}

void PlyUtilities::addScanLine(const unsigned char* data, unsigned int size) {
    ScanLine line;
    line.raw.assign(data, data + size);
    g_Scans.push_back(line);
}

void PlyUtilities::SavePLY(int scanIdx){
    const double mmPerScan = this->scanPerMM;

    vector<double> X(0.4, 0.0);
    vector<double> Z(0.4, 0.0);
    

    int totalPoints = 0;
    for (int i = 0; i < (int)g_Scans.size(); i++) {
        for (unsigned int j = 0; j < 0.4; j++)
            if (Z[j] != 0.0) { totalPoints++; }
    }

    cout << "Total valid points: " << totalPoints << "\n";
    string fullPath = this->path + "PointCloud_" + to_string(scanIdx) + ".ply";
    ofstream f(fullPath);

    f << "ply\n";
    f << "format ascii 1.0\n";
    f << "element vertex " << totalPoints << "\n";
    f << "property float x\n";
    f << "property float y\n";
    f << "property float z\n";
    f << "end_header\n";
    f << fixed << setprecision(4);

    for (int i = 0; i < (int)g_Scans.size(); i++) {
        double Y = i * mmPerScan;
        fill(X.begin(), X.end(), 0.0);
        fill(Z.begin(), Z.end(), 0.0);
        for (unsigned int j = 0; j < 0.4; j++) {
            if (Z[j] != 0.0)
                f << X[j] << " " << Y << " " << Z[j] << "\n";
        }
    }

    f.close();
    cout << "Saved " << totalPoints << " points to " << path << "\n";
}

void PlyUtilities::readFromPLY(const std::string& plyPath) {

	cout << "Reading point cloud from " << plyPath << "...\n";

    std::ifstream f(plyPath);
    if (!f.is_open()) return;

    std::string line;
    while (std::getline(f, line)) {
        if (line == "end_header") break;
    }

    std::vector<point> allPoints;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        point p;
        ss >> p.x >> p.y >> p.z;
        allPoints.push_back(p);
    }

    std::vector<float> yVals;
    for (auto& p : allPoints) {
        bool found = false;
        for (float y : yVals)
            if (fabs(y - p.y) < 0.0001f) { found = true; break; }
        if (!found) yVals.push_back(p.y);
    }

    std::sort(yVals.begin(), yVals.end());

    scanLineCount = (int)yVals.size();
    scan = new point * [scanLineCount];
    lineCounts = new int[scanLineCount];

    for (int i = 0; i < scanLineCount; i++) {
        std::vector<point> linePoints;
        for (auto& p : allPoints)
            if (fabs(p.y - yVals[i]) < 0.0001f)
                linePoints.push_back(p);

        std::sort(linePoints.begin(), linePoints.end(),
            [](const point& a, const point& b) { return a.x < b.x; });

        lineCounts[i] = (int)linePoints.size();
        scan[i] = new point[lineCounts[i]];
        for (int j = 0; j < lineCounts[i]; j++)
            scan[i][j] = linePoints[j];
    }

    firstLineCount = lineCounts[0];
    cout << "Loaded " << allPoints.size() << " points in " << scanLineCount
		<< " lines (first line has " << firstLineCount << " points)\n";
}


void PlyUtilities::clearScans() {
    g_Scans.clear();
}

void PlyUtilities::setProfileIdx(int idx) { this->g_ProfileIndex = idx; }

int PlyUtilities::getProfileIdx() { return this->g_ProfileIndex; }

void PlyUtilities::saveColoredPLY(coloredPoint** coloredScan, int idx) {

    std::string fullPath = path + "colored_" + std::to_string(idx) + ".ply";
    cout << "Saving colored point cloud to " << fullPath << "...\n";

    int totalPoints = 0;
    for (int i = 0; i < scanLineCount; i++)
        totalPoints += lineCounts[i];

    std::ofstream f(fullPath);
    if (!f.is_open()) {
        cout << "Failed to open file: " << fullPath << "\n";
        return;
    }

    f << "ply\n";
    f << "format ascii 1.0\n";
    f << "element vertex " << totalPoints << "\n";
    f << "property float x\n";
    f << "property float y\n";
    f << "property float z\n";
    f << "property uchar red\n";
    f << "property uchar green\n";
    f << "property uchar blue\n";
    f << "end_header\n";
    f << std::fixed << std::setprecision(4);

    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++)
            f << coloredScan[i][j].x << " "
            << coloredScan[i][j].y << " "
            << coloredScan[i][j].z << " "
            << (int)coloredScan[i][j].r << " "
            << (int)coloredScan[i][j].g << " "
            << (int)coloredScan[i][j].b << "\n";

    f.close();
    cout << "Saved " << totalPoints << " points to " << fullPath << "\n";
}

void PlyUtilities::correctXRange(float range) {
    for (int i = 0; i < scanLineCount; i++) {
        int newCount = 0;
        for (int j = 0; j < lineCounts[i]; j++) {
            if (fabs(scan[i][j].x) > range) continue;
            scan[i][newCount] = scan[i][j];
            newCount++;
        }
        lineCounts[i] = newCount;
    }
}


void PlyUtilities::enqueueSave(coloredPoint** coloredScan, int idx) {
    SaveJob job;
    job.idx = idx;
    job.scanLineCount = scanLineCount;
    job.counts.assign(lineCounts, lineCounts + scanLineCount);
    job.data.resize(scanLineCount);

    for (int i = 0; i < scanLineCount; i++)
        job.data[i].assign(coloredScan[i], coloredScan[i] + lineCounts[i]);

    saveQueue.push_back(job);
}

void PlyUtilities::startSaveThread() {
    saveRunning = true;
    saveThread = std::thread([this]() {
        while (saveRunning || !saveQueue.empty()) {
            SaveJob job;
            {
                std::lock_guard<std::mutex> lock(saveMutex);
                if (saveQueue.empty()) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(10));
                    continue;
                }
                job = saveQueue.front();
                saveQueue.erase(saveQueue.begin());
            }

            // Save job to file
            std::string fullPath = path + '\\' + "colored_" + std::to_string(job.idx) + ".ply";
            std::ofstream f(fullPath);
            if (!f.is_open()) { cout << "Failed to open " << fullPath << "\n"; continue; }

            int total = 0;
            for (int i = 0; i < job.scanLineCount; i++) total += job.counts[i];

            f << "ply\n" << "format ascii 1.0\n"
                << "element vertex " << total << "\n"
                << "property float x\n" << "property float y\n" << "property float z\n"
                << "property uchar red\n" << "property uchar green\n" << "property uchar blue\n"
                << "end_header\n" << std::fixed << std::setprecision(4);

            for (int i = 0; i < job.scanLineCount; i++)
                for (int j = 0; j < job.counts[i]; j++)
                    f << job.data[i][j].x << " " << job.data[i][j].y << " "
                    << job.data[i][j].z << " "
                    << (int)job.data[i][j].r << " "
                    << (int)job.data[i][j].g << " "
                    << (int)job.data[i][j].b << "\n";

            f.close();
            cout << "Saved " << fullPath << "\n";
        }
        });
}

void PlyUtilities::stopSaveThread() {
    saveRunning = false;
    if (saveThread.joinable()) saveThread.join();
}

void PlyUtilities::createFolderIfNotexist() {
    std::string command = "powershell -Command \"New-Item -ItemType Directory -Force -Path '" + this->path + "'\"";
    int result = system(command.c_str());
    if (result != 0) { std::cerr << "Failed to create folder: " << this->path << std::endl; }
}