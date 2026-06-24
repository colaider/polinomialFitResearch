#include "HyperME.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;
using point = PlyUtilities::point;
using coloredPoint = PlyUtilities::coloredPoint;

vector<point> HyperME::findRefernceLine(float region) {
    float minZ = FLT_MAX, maxZ = -FLT_MAX;

    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++) {
            if (scan[i][j].z < minZ) minZ = scan[i][j].z;
            if (scan[i][j].z > maxZ) maxZ = scan[i][j].z;
        }

    float threshold = minZ + (maxZ - minZ) * region;

    vector<point> pts;
    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++)
            if (scan[i][j].z >= threshold)
                pts.push_back(scan[i][j]);

    cout << "findReferenceLine: kept " << pts.size() << " points above " << threshold << "\n";
    return pts;
}

void HyperME::fitHyperSurface(vector<point> pts, int degree) {

    const int N = 15;
    double M[N][N] = {};
    double b[N] = {};

    normCx = 0; normCy = 0; normCz = 0;
    normSx = 1; normSy = 1; normSz = 1;

    for (int i = 0; i < N; i++) coeffs[i] = 0;

    // Define which coefficients belong to which degree
    // phi index: 0-4=degree4, 5-8=degree3, 9-11=degree2, 12-13=degree1, 14=degree0
    int coeffDegree[N] = { 4,4,4,4,4, 3,3,3,3, 2,2,2, 1,1, 0 };

    // Build mask - which coefficients to use
    bool use[N];
    for (int i = 0; i < N; i++)
        use[i] = (coeffDegree[i] <= degree);

    for (int i = 0; i < (int)pts.size(); i++) {
        normCx += pts[i].x;
        normCy += pts[i].y;
        normCz += pts[i].z;
    }
    normCx /= pts.size(); normCy /= pts.size(); normCz /= pts.size();

    float maxXs = -FLT_MAX, maxYs = -FLT_MAX, maxZs = -FLT_MAX;
    for (int i = 0; i < (int)pts.size(); i++) {
        if (maxXs < fabs(pts[i].x - normCx)) maxXs = fabs(pts[i].x - normCx);
        if (maxYs < fabs(pts[i].y - normCy)) maxYs = fabs(pts[i].y - normCy);
        if (maxZs < fabs(pts[i].z - normCz)) maxZs = fabs(pts[i].z - normCz);
    }
    normSx = maxXs; normSy = maxYs; normSz = maxZs;
    if (normSx < 1e-6f) normSx = 1.0f;
    if (normSy < 1e-6f) normSy = 1.0f;
    if (normSz < 1e-6f) normSz = 1.0f;

    for (int i = 0; i < (int)pts.size(); i++) {
        double x = (pts[i].x - normCx) / normSx;
        double y = (pts[i].y - normCy) / normSy;
        double z = (pts[i].z - normCz) / normSz;

        double phi[N] = {
            x * x * x * x, y * y * y * y, x * x * x * y, y * y * y * x, x * x * y * y,
            x * x * x, y * y * y, x * x * y, x * y * y,
            x * x, y * y, x * y,
            x, y, 1.0
        };

        // zero out unused terms
        for (int k = 0; k < N; k++)
            if (!use[k]) phi[k] = 0.0;

        for (int r = 0; r < N; r++) {
            if (!use[r]) continue;
            b[r] += phi[r] * z;
            for (int c = r; c < N; c++) {
                if (!use[c]) continue;
                M[r][c] += phi[r] * phi[c];
            }
        }
    }

    // Set diagonal to 1 for unused coefficients to avoid singular matrix
    for (int i = 0; i < N; i++)
        if (!use[i]) M[i][i] = 1.0;

    for (int r = 0; r < N; r++)
        for (int c = 0; c < r; c++)
            M[r][c] = M[c][r];

    for (int col = 0; col < N; col++) {
        int maxRow = col;
        for (int row = col + 1; row < N; row++)
            if (fabs(M[row][col]) > fabs(M[maxRow][col])) maxRow = row;
        for (int k = 0; k < N; k++) swap(M[col][k], M[maxRow][k]);
        swap(b[col], b[maxRow]);

        for (int row = col + 1; row < N; row++) {
            double f = M[row][col] / M[col][col];
            for (int k = col; k < N; k++)
                M[row][k] -= f * M[col][k];
            b[row] -= f * b[col];
        }
    }

    for (int row = N - 1; row >= 0; row--) {
        double sum = b[row];
        for (int k = row + 1; k < N; k++)
            sum -= M[row][k] * coeffs[k];
        coeffs[row] = use[row] ? (float)(sum / M[row][row]) : 0.0f;
    }

    std::cout << "Degree " << degree << " Coefficients: ";
    for (int i = 0; i < N; i++) std::cout << coeffs[i] << " ";
    std::cout << "\n";
}

float HyperME::eval(float x, float y) {
    float xn = (x - normCx) / normSx;
    float yn = (y - normCy) / normSy;

    float zn = coeffs[0] * xn * xn * xn * xn + coeffs[1] * yn * yn * yn * yn +
        coeffs[2] * xn * xn * xn * yn + coeffs[3] * xn * yn * yn * yn +
        coeffs[4] * xn * xn * yn * yn +
        coeffs[5] * xn * xn * xn + coeffs[6] * yn * yn * yn +
        coeffs[7] * xn * xn * yn + coeffs[8] * xn * yn * yn +
        coeffs[9] * xn * xn + coeffs[10] * yn * yn +
        coeffs[11] * xn * yn + coeffs[12] * xn +
        coeffs[13] * yn + coeffs[14];

    return zn * normSz + normCz;
}

vector<point> HyperME::getSubCloud(float zMin, float zMax) {
	vector<point> subCloud;
    float maxY = -FLT_MAX;
    float minY = FLT_MAX;
    float maxZ = -FLT_MAX;
    float minZ = FLT_MAX;
    float treshold = 0.0f;

    point maxPoint, minPoint, rMaxPoint, rMinPoint;

    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++) {
            if (scan[i][j].z < minZ) minZ = scan[i][j].z;
            if (scan[i][j].z > maxZ) maxZ = scan[i][j].z;
        }


    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++) {
            if (scan[i][j].z <= treshold) continue;
            if (scan[i][j].y < minY) minY = scan[i][j].y;
            if (scan[i][j].y > maxY) maxY = scan[i][j].y;
        }


    for (int i = 0; i < scanLineCount; i++) {
        for (int j = 0; j < lineCounts[i]; j++) {
            float x = scan[i][j].x;
            float y = scan[i][j].y;
            float z = scan[i][j].z;
            float zSurface = eval(x, y);
            float dist = z - zSurface;
			float yHat = (y - minY) / (maxY - minY);
			float group = zMin + (zMax - zMin) * yHat;

            if (dist < -group) continue;
            subCloud.push_back(scan[i][j]);
        }
    }
	return subCloud;
}


void HyperME::colorByPlaneDistance() {
    float range = 0.0f;
	fitHyperSurface(findRefernceLine(range), 2);
	int degree = 2; float band = 0.2, bandMax = 0.32, bandMin = 0.05; float initialCurvature = -FLT_MAX;

    auto fit = [&](float& m, float& b, vector<float> x, vector<float> y) {
        float sumOfY = 0.0f, sumOfX = 0.0f, sumOfX2 = 0.0f, sumOfXY = 0.0f;
        int n = 0;
        for (int i = 0; i < x.size(); i++) {
            sumOfY += y[i]; sumOfX += x[i];
            sumOfX2 += x[i] * x[i]; sumOfXY += x[i] * y[i]; n++;
        }
        m = (n * sumOfXY - sumOfX * sumOfY) / (n * sumOfX2 - sumOfX * sumOfX);
        b = (sumOfY * sumOfX2 - sumOfX * sumOfXY) / (n * sumOfX2 - sumOfX * sumOfX);
        };

    auto curvToBand = [&](float x) {
        vector<float> X = { 0.0005, 0.005, 0.01, 0.1, 0.5 };
		vector<float> Y = { 0.05, 0.1, 0.15, 0.25, 0.3 };
        Y = this->hetIsHoles() ? vector<float> { 0.15, 0.18, 0.2, 0.25, 0.3 } : Y;
        for (size_t k = 0; k < X.size(); k++) { X[k] = log10(X[k]); }
		float m = 0, b = 0; fit(m, b, X, Y);

        std::cout << "m,n::" << m << "," << b << ":\n";
        
        float result = m * log10(max(x, 1e-6f)) + b; 
		result = result < bandMin ? bandMin : result;
		result = result > bandMax ? bandMax : result;
        return result;
    };    

	auto fitSurfceOfDegree = [&](int degree, float band, int ji, bool lastIter) {
	    ji = ji== -1 ? (int)(0.2 / band) : ji;
        fitHyperSurface(getSubCloud(0.4 / ji, 0.6 / ji ), degree);

        coloredScan = new PlyUtilities::coloredPoint * [scanLineCount];

        for (int i = 0; i < scanLineCount; i++) {
            coloredScan[i] = new PlyUtilities::coloredPoint[lineCounts[i]];

            for (int j = 0; j < lineCounts[i]; j++) {
                float x = scan[i][j].x;
                float y = scan[i][j].y;
                float z = scan[i][j].z;

                float zSurface = eval(x, y);

                float dist = fabs(z - zSurface);

                coloredScan[i][j].x = x;
                coloredScan[i][j].y = y;
                coloredScan[i][j].z = z;

                bool condition = lastIter? z >= zSurface || dist <= band : dist <= band;
                if (condition) {
                    coloredScan[i][j].r = 0;
                    coloredScan[i][j].g = 255;
                    coloredScan[i][j].b = 0;
                }
                else {
                    coloredScan[i][j].r = 255;
                    coloredScan[i][j].g = 0;
                    coloredScan[i][j].b = 0;
                }
            }
        }
		};

    auto rateOfCChange = [&](int d1, int d2, float band) {
        fitSurfceOfDegree(d1, band, -1, false);
        float secDeg = avgGaussianCurvature();
        fitHyperSurface(findRefernceLine(range), d1);
        fitSurfceOfDegree(d2, band, -1, false);
        float thirdDeg = avgGaussianCurvature();
        float dC = (100 * secDeg) / (100 * thirdDeg);
        fitHyperSurface(findRefernceLine(range), d2);
        return dC;
        };

    auto polynomialDegree = [&](float band) {
        float alpha = 0.6;
        float dC1 = rateOfCChange(2, 3, band);
        float dC2 = rateOfCChange(3, 4, band);

        int predictedDegree = dC1 <= alpha ? 2 : 3;
        if (predictedDegree == 3 && dC2 > alpha) predictedDegree = 4;
        cout << "predicted degree: " << predictedDegree << ", dC1: " << dC1 << ", dC2: " << dC2 << "\n";
        return predictedDegree;
        };


	degree = polynomialDegree(band);

    fitHyperSurface(findRefernceLine(range), degree);
    
    for (int i = 1; i <= 2; i++) {
        fitSurfceOfDegree(degree, band, i, i < 2 ? false : true);
        float curv = avgGaussianCurvature();
        band = 0.08f;
		std::cout << "band: " << band << "\n";
    }
}


void HyperME::flattenPointCloud() {
    for (int i = 0; i < scanLineCount; i++) {
        for (int j = 0; j < lineCounts[i]; j++) {
            float x = coloredScan[i][j].x;
            float y = coloredScan[i][j].y;

            float zSurface = eval(x, y);
			coloredScan[i][j].z = coloredScan[i][j].z - zSurface;
            
        }
    }
}

float HyperME::scanDisplacementCorrrection(float range) {

    float globalMinX = FLT_MAX;
    float lastLineMinX = FLT_MAX;
    float maxY = -FLT_MAX, minY = FLT_MAX;

    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++)
            if (coloredScan[i][j].g == 255 && coloredScan[i][j].r == 0) {
                if (coloredScan[i][j].x < globalMinX) globalMinX = coloredScan[i][j].x;
                if (coloredScan[i][j].y > maxY)       maxY = coloredScan[i][j].y;
                if (coloredScan[i][j].y < minY)       minY = coloredScan[i][j].y;
            }

    float targetY = maxY - blueRegion - 0.5f;

    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++)
            if (coloredScan[i][j].g == 255 && coloredScan[i][j].r == 0)
                if (fabs(coloredScan[i][j].y - targetY) < 0.5f)
                    if (coloredScan[i][j].x < lastLineMinX) lastLineMinX = coloredScan[i][j].x;

    std::cout << "globalMinX=" << globalMinX << " lastLineMinX=" << lastLineMinX << "\n";

    if (globalMinX == FLT_MAX || lastLineMinX == FLT_MAX) return 0.0f;
    float diff = lastLineMinX - globalMinX;

    for (int i = 0; i < scanLineCount;i++)
        for (int j = 0;j < lineCounts[i];j++)
            coloredScan[i][j].y -= minY;

    if (fabs(diff) >= 0.05f) return fabs(diff) * cos(rotTeta);
    return 0.0f;
}

float HyperME::scanDisplacementCorrrectionHoles() {
    // Step 1 - global rightmost green point
    float globalMinX = FLT_MAX;
    for (int i = 0; i < scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++)
            if (coloredScan[i][j].g == 255 && coloredScan[i][j].r == 0)
                if (coloredScan[i][j].x < globalMinX) globalMinX = coloredScan[i][j].x;

    // Find middle 60% scan line range
    int startLine = (int)(scanLineCount * 0.2f);
    int endLine = (int)(scanLineCount * 0.8f);

    // Step 2 - per scan line smallest X green point in middle 60%
    // Step 3 - find furthest from globalMinX
    float maxDiff = 0.0f;

    for (int i = startLine; i < endLine; i++) {
        float lineMinX = FLT_MAX;
        for (int j = 0; j < lineCounts[i]; j++)
            if (coloredScan[i][j].g == 255 && coloredScan[i][j].r == 0)
                if (coloredScan[i][j].x < lineMinX) lineMinX = coloredScan[i][j].x;

        if (lineMinX == FLT_MAX) continue;

        float diff = fabs(lineMinX - globalMinX);
        if (diff > maxDiff) maxDiff = diff;
    }

    cout << "globalMinX=" << globalMinX << " maxDiff=" << maxDiff
        << " startLine=" << startLine << " endLine=" << endLine << "\n";

    if (maxDiff < 0.05f) return 0.0f;
    return maxDiff * cos(rotTeta);
}

void HyperME::alignYAxis() {

    float maxY = -FLT_MAX, minY = FLT_MAX;

    for(int i=0; i<scanLineCount; i++)
        for (int j = 0; j < lineCounts[i]; j++) {
            if (coloredScan[i][j].g == 0 || coloredScan[i][j].r == 255) continue;
            if (coloredScan[i][j].y > maxY && coloredScan[i][j].g) maxY = coloredScan[i][j].y;
            if (coloredScan[i][j].y < minY) minY = coloredScan[i][j].y;
        }

    float bottomTreshold = minY + (maxY - minY) * 0.15;
	float topTreshold = maxY - (maxY - minY) * 0.15;

    vector<coloredPoint> rightPoints;
     
    for (int i = 0; i < scanLineCount; i++) {
		coloredPoint maxX; maxX.x = -FLT_MAX;
        for(int j=0; j<lineCounts[i]; j++) 
            if (coloredScan[i][j].x > maxX.x && coloredScan[i][j].g == 255 && coloredScan[i][j].r == 0){
                if(coloredScan[i][j].y < topTreshold && coloredScan[i][j].y > bottomTreshold)
                    maxX = coloredScan[i][j];
			}
        if (maxX.x == -FLT_MAX) continue;
		//cout << "Line " << i << ": maxX = (" << maxX.x << ", " << maxX.y << ", " << maxX.z << ")\n";
		rightPoints.push_back(maxX);
    }


    if (rightPoints.size() < 4) {
		cout << "Not enough points to align Y axis\n"; return;
    }

	int quter = (int)(rightPoints.size() / 4);
    vector<coloredPoint> filterd;

    for (int i = 0; i < quter; i++){
        float max = -FLT_MAX; int maxJ = 0;
        for (int j = 0; j <  rightPoints.size(); j++) 
			if (rightPoints[j].x > max) { max = rightPoints[j].x; maxJ = j; }
            
        filterd.push_back(rightPoints[maxJ]);
        rightPoints.erase(rightPoints.begin() + maxJ);
    }



    float sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    int   n = (int)filterd.size();

    for (auto& p : filterd) {
        sumX += p.x;
        sumY += p.y;
        sumXY += p.x * p.y;
        sumX2 += p.x * p.x;
    }

    float m = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    float b = (sumY - m * sumX) / n;

	cout << "Reference line: y = " << m << "*x + " << b << "\n";

    float angle = rotTeta = (PI/2.0f)+atan(m);
    float cosA = cos(-angle);
    float sinA = sin(-angle);
    dYScaling = cosA;
    float denom = sqrt(m * m + 1.0f);

    cout << "Rotation angle: " << angle * 180.0f / 3.14159f << " degrees\n";


    for (int i = 0; i < scanLineCount; i++) {
        for (int j = 0; j < lineCounts[i]; j++) {
            float x = coloredScan[i][j].x;
            float y = coloredScan[i][j].y;

            float xNew = x * cosA - y * sinA;
            float yNew = x * sinA + y * cosA;

            coloredScan[i][j].x = xNew;  // modify in place
            coloredScan[i][j].y = yNew;

            if (fabs(y - m * x - b) / denom <= 0.2f) {
                coloredScan[i][j].r = 255;
                coloredScan[i][j].g = 105;
                coloredScan[i][j].b = 180;
            }
        }
    }
    cout << "Rotated and clipped point cloud\n";

    plyUtils->setCounts(lineCounts, scanLineCount);
}

void HyperME::applyRotation(float teta) {
    float cosA = cos(-teta);
    float sinA = sin(-teta);
    dYScaling = cosA;
	rotTeta = teta;

    cout << "Rotation angle: " << teta * 180.0f / 3.14159f << " degrees\n";

    for (int i = 0; i < scanLineCount; i++) {
        for (int j = 0; j < lineCounts[i]; j++) {
            float x = coloredScan[i][j].x;
            float y = coloredScan[i][j].y;

            float xNew = x * cosA - y * sinA;
            float yNew = x * sinA + y * cosA;

            coloredScan[i][j].x = xNew;  
            coloredScan[i][j].y = yNew;
        }
    }
    cout << "Rotated and clipped point cloud\n";

    plyUtils->setCounts(lineCounts, scanLineCount);
}

void HyperME::alignByEigenVector() {
    float xC = 0.0f, yC = 0.0f; int n = 0;
    float sXsX = 0.0f, sYsY = 0.0f, sXsY = 0.0f;

    for (int i = 0;i < scanLineCount;i++) 
        for (int j = 0;j < lineCounts[i];j++) { 
            if (coloredScan[i][j].x < 10.0f ||(coloredScan[i][j].g == 0)) continue;
            xC += coloredScan[i][j].x; 
            yC += coloredScan[i][j].y;n++; 
        }

    xC /= n; yC /= n;


    for (int i = 0;i < scanLineCount;i++)
        for (int j = 0;j < lineCounts[i];j++){
            if (coloredScan[i][j].x < 10.0f || (coloredScan[i][j].g == 0)) continue;
			sXsX += (coloredScan[i][j].x - xC) * (coloredScan[i][j].x - xC);
            sYsY += (coloredScan[i][j].y - yC) * (coloredScan[i][j].y - yC);
            sXsY += (coloredScan[i][j].y - yC) * (coloredScan[i][j].x - xC);
         }
    sXsX /= n; sYsY /= n; sXsY /= n;

    vector<vector<float>> C;
    vector<float> cRow1; vector<float> cRow2;
	cRow1.push_back(sXsX); cRow1.push_back(sXsY); C.push_back(cRow1);
	cRow2.push_back(sXsY); cRow2.push_back(sYsY); C.push_back(cRow2);
    float a = 1.0f, b = -(sXsX + sYsY), c = sXsX * sYsY - sXsY * sXsY;

    float D = b * b - 4 * a * c;
    float l1 = (-b + sqrt(D)) / (2 * a), l2 = (-b - sqrt(D)) / (2 * a);

    float m = (l1 < l2) ? (sXsY / (l1 - sXsX)) : (sXsY / (l2 - sXsX)); b = 0;


    float angle = rotTeta = atan(1.0f / m);
    float cosA = cos(-angle);
    float sinA = sin(-angle);
    dYScaling = cosA;
    float denom = sqrt(m * m + 1.0f);

    cout << "Rotation angle: " << angle * 180.0f / 3.14159f << " degrees\n";

    for (int i = 0; i < scanLineCount; i++) {
        for (int j = 0; j < lineCounts[i]; j++) {
            float x = coloredScan[i][j].x;
            float y = coloredScan[i][j].y;

            float xNew = x * cosA - y * sinA;
            float yNew = x * sinA + y * cosA;

            coloredScan[i][j].x = xNew;  // modify in place
            coloredScan[i][j].y = yNew;

            if (fabs(y - m * x - b) / denom <= 0.2f) {
                coloredScan[i][j].r = 255;
                coloredScan[i][j].g = 105;
                coloredScan[i][j].b = 180;
            }
        }
    }
    cout << "Rotated and clipped point cloud\n";

    plyUtils->setCounts(lineCounts, scanLineCount);
}

void HyperME::alignXAxis() {
    float miX = FLT_MAX; float minY0 = FLT_MAX; float minYX = FLT_MAX;

    for (int i = 0; i < scanLineCount;i++)
        for (int j = 0;j < lineCounts[i];j++) {
            if (coloredScan[i][j].g != 255 || coloredScan[i][j].r != 0) continue;
            if (coloredScan[i][j].x < miX) miX = coloredScan[i][j].x;
        }
            
    for (int i = 0;i < scanLineCount;i++) {
        for (int j = 0;j < lineCounts[i]; j++) {
            if (coloredScan[i][j].g != 255 || coloredScan[i][j].r != 0) continue;
            if(fabs(coloredScan[i][j].x - miX) < 0.2f)
                if(coloredScan[i][j].y < minYX) minYX = coloredScan[i][j].y;
            if(fabs(coloredScan[i][j].x) < 0.2f)
                if(coloredScan[i][j].y < minY0) minY0 = coloredScan[i][j].y;    
        }
    }

    auto fit = [&](float& m, float& b,  float r) {
        float sumOfY = 0.0f, sumOfX = 0.0f, sumOfX2 = 0.0f, sumOfXY = 0.0f;
        int n = 0;
        for (int i = 0; i < scanLineCount;i++) {
            for (int j = 0; j < lineCounts[i]; j++) {
                if (coloredScan[i][j].g != 255 || coloredScan[i][j].r != 0) continue;
                if (coloredScan[i][j].y <= (coloredScan[i][j].x * m + b) + r) {
                    sumOfY += coloredScan[i][j].y; sumOfX += coloredScan[i][j].x;
                    sumOfX2 += coloredScan[i][j].x * coloredScan[i][j].x;
                    sumOfXY += coloredScan[i][j].x * coloredScan[i][j].y;
                    n++;
                }
            }
        }
        m = (n * sumOfXY - sumOfX * sumOfY) / (n * sumOfX2 - sumOfX * sumOfX);
        b = (sumOfY * sumOfX2 - sumOfX * sumOfXY) / (n * sumOfX2 - sumOfX * sumOfX);

        };

    float m = (minY0 - minYX) / (0 - miX); float b = minY0;

    for (float r : {1.0f, 0.8f, 0.5f, 0.1f, 0.08f, 0.06f, 0.04f, 0.02f}) fit(m, b, r);


    float angle = rotTeta = atan(m);
    float cosA = cos(-angle);
    float sinA = sin(-angle);
    dYScaling = cosA;

    std::cout << "Rotation angle: " << angle * 180.0f / 3.14159f << " degrees\n";


    for (int i = 0; i < scanLineCount; i++) {
        for (int j = 0; j < lineCounts[i]; j++) {
            float x = coloredScan[i][j].x;
            float y = coloredScan[i][j].y;

            float xNew = x * cosA - y * sinA;
            float yNew = x * sinA + y * cosA;

            coloredScan[i][j].x = xNew;  // modify in place
            coloredScan[i][j].y = yNew;

            if (fabs(y - (m * x + b)) <= 0.3f) {
                coloredScan[i][j].r = 255;
                coloredScan[i][j].g = 105;
                coloredScan[i][j].b = 180;
            }
        }
    }
    std::cout << "Rotated and clipped point cloud\n";

    plyUtils->setCounts(lineCounts, scanLineCount);
}

float HyperME::gaussianCurvature(float xRaw, float yRaw) {
	float x = (xRaw - normCx) / normSx; float y = (yRaw - normCy) / normSy;

	float f_x = 4*coeffs[0]*x*x*x + 3*coeffs[2]*x*x*y + coeffs[3]*y*y*y + 2*coeffs[4]*y*y*x + 3*coeffs[5]*x*x + 
                2*coeffs[7]*x*y + coeffs[8]*y*y + 2*coeffs[9]*x + coeffs[11]*y + coeffs[12];
    float f_y = 4*coeffs[1]*y*y*y + 3*coeffs[3]*x*y*y + coeffs[2]*x*x*x + 2*coeffs[4]*x*x*y + 3*coeffs[6]*y*y + 
                2*coeffs[8]*x*y + coeffs[7]*x*x + 2*coeffs[10]*y + coeffs[11]*x + coeffs[13];

    float f_xx = 12*coeffs[0]*x*x + 6*coeffs[2]*x*y + 2*coeffs[4]*y*y + 6*coeffs[5]*x + 2*coeffs[7]*y + 2*coeffs[9];
    float f_yy = 12*coeffs[1]*y*y + 6*coeffs[3]*x*y + 2*coeffs[4]*x*x + 6*coeffs[6]*y + 2*coeffs[8]*x + 2*coeffs[10];
    float f_xy = 3*coeffs[2]*x*x + 3*coeffs[3]*y*y + 4 *coeffs[4]*x*y + 2*coeffs[7]*x + 2*coeffs[8]*y + coeffs[11];

    float firstFormDet = (f_x * f_x + f_y * f_y + 1) * (f_x * f_x + f_y * f_y + 1);
	float secondFormDet = (f_xx * f_yy) - (f_xy * f_xy);

    if(fabs(firstFormDet) < 1e-6) return 0.0f;
    float K = secondFormDet/firstFormDet;
    return K;
}

float HyperME::avgGaussianCurvature() {
	float avgK = 0.0f; int count = 0;
    for(int i = 0; i < scanLineCount; i++) {
        for(int j = 0; j < lineCounts[i]; j++) {
            if(coloredScan[i][j].g != 255 || coloredScan[i][j].r != 0) continue;
            avgK += fabs(gaussianCurvature(coloredScan[i][j].x, coloredScan[i][j].y));
            count++;
        }
    }
	if (count == 0) return 0.0f;
    std::cout << "\n" << "avg Curvature =" << avgK / count << "\n";
    return avgK / count;
}