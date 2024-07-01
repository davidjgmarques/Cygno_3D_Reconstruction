#include <TApplication.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TLine.h>
#include <TRandom3.h>

int main(int argc, char** argv) {
    // Create a ROOT application
    TApplication app("example", &argc, argv);

    // Create a canvas to draw on
    TCanvas canvas("canvas", "3D Points Example", 800, 600);

    // Define the dimensions of the box
    const Int_t nBins = 30;
    const Double_t xMin = 0, xMax = 33;
    const Double_t zMin = 0, zMax = 33;
    const Double_t yMin = 0, yMax = 50;

    // Create a TH3F object to define the 3D box and hold the points
    TH3F *h3 = new TH3F("h3", "3D Points", 33, xMin, xMax, 50, yMin, yMax, 33, zMin, zMax);

    // Define parameters for the virtual line
    Double_t lineX1 = 10, lineY1 = 30, lineZ1 = 10; // Starting point of the line
    Double_t lineX2 = 10, lineY2 = 20, lineZ2 = 25;   // Ending point of the line

    // Generate random points around the virtual line
    TRandom3 randGen; // Random number generator
    for (Int_t i = 0; i < 100; ++i) {
        // Generate random point coordinates around the line
        Double_t t = randGen.Uniform(0.0, 1.0); // Parameter along the line
        Double_t x = lineX1 + t * (lineX2 - lineX1) + randGen.Gaus(0.0, 0.06); // Add Gaussian noise
        Double_t y = lineY1 + t * (lineY2 - lineY1) + randGen.Gaus(0.0, 0.06); // Add Gaussian noise
        Double_t z = lineZ1 + t * (lineZ2 - lineZ1) + randGen.Gaus(0.0, 0.06); // Add Gaussian noise

        // Fill the point into the histogram
        h3->Fill(x, y, z);
    }

    // Draw the 3D histogram with grid lines
    h3->SetMarkerColor(kOrange+10);
    h3->SetMarkerStyle(20);
    h3->Draw("");  // "BOX" option draws the 3D histogram with axis grids

    // Draw the virtual line
    TLine *line = new TLine(lineX1, lineY1, lineX2, lineY2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->Draw("same");

    // Set the range of the axes
    h3->GetXaxis()->SetRangeUser(xMin, xMax);
    h3->GetYaxis()->SetRangeUser(yMin, yMax);
    h3->GetZaxis()->SetRangeUser(zMin, zMax);

    h3->GetXaxis()->SetTitle("X");
    h3->GetYaxis()->SetTitle("Z");
    h3->GetZaxis()->SetTitle("Y");

    // Set colors for each axis
    h3->GetXaxis()->SetAxisColor(kRed);
    h3->GetYaxis()->SetAxisColor(kGreen+1);
    h3->GetZaxis()->SetAxisColor(kBlue+1);
    h3->GetXaxis()->SetLabelColor(kRed);
    h3->GetYaxis()->SetLabelColor(kGreen+1);
    h3->GetZaxis()->SetLabelColor(kBlue+1);

    // Update the canvas
    // canvas.Update();

    // Run the application
    app.Run();

    return 0;
}
