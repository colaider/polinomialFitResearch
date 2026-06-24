#define _CRT_SECURE_NO_WARNINGS
#define NOMINMAX
#include <iostream>
#include <string>
#include "HyperME.h"
#include "PlyUtilities.h"
#include "MyMetrologyEngine.h"

using namespace std;
int main(int argc, char* argv[]) {
  
    // ---------------------------------------------------------------
    // 1. Configuration — edit these before running
    // ---------------------------------------------------------------
    string plyPath = "./ply/output.ply"; 
    string sampleName = "test_sample";
    string outputFolder = "./output";

    float sectionY = 50.0f;   
    float sectionX = 20.0f;   

    PlyUtilities plyUtils(outputFolder, sampleName, 40);
    cout << "Loading PLY from: " << plyPath << "\n";
    plyUtils.readFromPLY(plyPath);
    
    if (plyUtils.getScanLineCount() == 0) {
        cerr << "ERROR: No scan lines loaded. Check file path.\n";
        return -1;
    }
    cout << "Loaded " << plyUtils.getScanLineCount() << " scan lines.\n";

    // ---------------------------------------------------------------
    // 3. Run the metrology pipeline
    // ---------------------------------------------------------------
    HyperME mEngine(plyPath, &plyUtils);
    mEngine.setRange(sectionY);
  
    mEngine.getPlyUtils()->correctXRange(sectionX / 2.0f);
    // Fit polynomial surface and colour green/red by distance
    cout << "Running colorByPlaneDistance...\n";
    mEngine.colorByPlaneDistance();
    // pause here first
    // ---------------------------------------------------------------
    // 5. Save colored PLY
    // ---------------------------------------------------------------
    plyUtils.startSaveThread();
    plyUtils.enqueueSave(mEngine.getCloloredScan(), 0);
    plyUtils.stopSaveThread();
    cout << "Colored PLY saved to: " << outputFolder << "\\" << sampleName << "\n";

    cout << "\nDone. Press Enter to exit.\n";
    cin.get();
    return 0;
}