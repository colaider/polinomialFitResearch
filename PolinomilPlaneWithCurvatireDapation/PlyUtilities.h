#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <chrono>
#include <winsock2.h>

class PlyUtilities {
public:
    struct point {
        float x, y, z;
    };

    struct coloredPoint : public point {
        unsigned char r, g, b;
    };;

private:
    std::string          path;
    double               scanPerMM;

    struct ScanLine { std::vector<unsigned char> raw; };
    std::vector<ScanLine> g_Scans;
    int g_ProfileIndex = 0;

    point**  scan = nullptr;
    int*     lineCounts = nullptr;
    int      scanLineCount = 0;
    int      firstLineCount = 0;
    int      lastScanLineCount = 0;
	std::string  sampleName;

    std::mutex            saveMutex;
    std::thread           saveThread;
    bool                  saveRunning = false;

public:

    struct SaveJob {
        std::vector<std::vector<coloredPoint>> data;  
        std::vector<int> counts;
        int scanLineCount;
        int idx;
    };

    std::vector<SaveJob> saveQueue;

    PlyUtilities(std::string path, std::string sampleName, int devisor = 40);
    void addScanLine(const unsigned char* data, unsigned int size);
    void SavePLY(int scanIdx);
    void clearScans();
    void setProfileIdx(int idx);
    int  getProfileIdx();
    void readFromPLY(const std::string& plyPath);
    void readColoredFromPLY(const std::string& plyPath);
    void saveColoredPLY(coloredPoint** coloredScan, int idx);
    void writeTextToCloud(std::string text, float x, float y);
    void correctXRange(float range);
    void enqueueSave(coloredPoint** coloredScan, int idx);
    void startSaveThread();
    void stopSaveThread();
    void saveFinalPLY(const std::string& sourceFolder, const std::string& versionName, const std::string& outFolder);
    void createFolderIfNotexist();


    point** getScan() { return scan; }
    int*    getLineCounts() { return lineCounts; }
    int     getScanLineCount() { return scanLineCount; }
    int     getFirstLineCount() { return firstLineCount; }
    void    setCounts(int* lc, int slc) { lineCounts = lc; scanLineCount = slc; }
    double  getScanPerMM() { return scanPerMM; };
	void    setDensity(float density) { scanPerMM = density / 1000.0; }

};