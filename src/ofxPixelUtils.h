#pragma once

#include "ofMain.h"

// begin namespace
namespace ofxPixelUtils {


/// Color functions:

float getMeanColor(const ofPixels& pixels, int channel);

vector<float> getMeanColor(const ofPixels& pixels);

float getBrightness(const ofPixels& pixels);

void adjustGain(ofPixels& pixels, float gain);
void adjustGain(ofPixels& pixels, float gain, int channel);

// Adjust brightness and contrast.
// Contrast is a factor, brightness is an offset.
// Setting the contrast to 1 will leave the contrast unchanged (uses no CPU).
// Same is true for a brightness of 0
void adjustBrightnessAndContrast(ofPixels& pixels, float brightness, float contrast);
void adjustBrightnessAndContrast(ofPixels& pixels, float brightness, float contrast, int channel);

// Adjust hue, saturation and brightness with float factors (only for RGB and RGBA pixels)
// hue is an offset where -0.5f and 0.5f correspond to -180 degree and +180 degree in the color circle
// saturation and brightness are factors
void adjustHSB(ofPixels& pixels, float hue, float saturation, float brightness);

// Apply color curves (represented as unsigned char vectors).
void applyColorCurve(ofPixels& pixels, vector<unsigned char>& table, int channel);
void applyColorCurve(ofPixels& pixels, vector<unsigned char>& table);

// convert from RBG to greyscale
enum class GrayscaleMode{
	LIGHTNESS, 
	AVERAGE, 
	LUMINANCE
};

void convertToGrayscale(const ofPixels& src, ofPixels& dest, GrayscaleMode mode = GrayscaleMode::AVERAGE);

// void invertColor(ofPixels& pixels, int channel);
// void invertColor(ofPixels& pixels);

// void applyThreshold(ofPixels& pixels, float low, float high, int channel, bool invert); // pass color only if it's within the range of [low, high] (invert == true -> opposite behavior)

// void applyAlphaMask(ofPixels& pixels, float redLow, float redHigh, float greenLow, float greenHigh, float blueLow, float blueHigh, bool mode);

// ofPixels convolveWithKernel(const ofPixels& pix, const vector<float>& kernel);

void applyAlphaMask(const ofPixels& src, const ofPixels& mask, ofPixels& dest, bool invertMask = false);

/// Movement detection

class MovementDetection {
protected:
    float *buf_1, *buf_w;
    float myFeedback, myGain;
    int myThreshold;
    bool bBinary;
    bool bInvert;
    ofPixels outPixels;
    int width;
    int height;
    int numChannels;
    bool bAllocated;
public:
    /// constructors
    MovementDetection();
    MovementDetection(int w, int h, int channels);
    /// destructor
    ~MovementDetection();
    /// allocate (will be called automatically by MovementDetection::in
    /// at the first time or when dimensions change
    void allocate(int w, int h, int channels);
    /// ask if buffers are allocated
    bool isAllocated();
    /// send ofPixels and calculate new filter state
    void in(const ofPixels& inPixels);
    /// get current filter state as ofPixels reference
    const ofPixels& out() const;
    /// clear the filter
    void clearFilter();

    void setAveraging(float coeff);
    void setGain(float gain);
    void setBinary(bool mode);
    void setThreshold(float threshold);
    void setInvert(bool mode);
    float getAveraging() const;
    float getGain() const;
    bool isBinary() const;
    float getThreshold() const;
    bool isInvert() const;
};

/// Background subtraction

class BackgroundSubstraction {
protected:
    ofPixels myBackground;
    ofPixels outPixels;
    int myThreshold;
    bool bBinary;
    bool bInvert;
    int width;
    int height;
    int numChannels;
public:
    BackgroundSubstraction();
    ~BackgroundSubstraction();
    void setBackground(const ofPixels& pix);
    void setBackground(int w, int h, int channels, ofColor color);
    const ofPixels& getBackground() const;
    bool isBackgroundSet() const;
    void in(const ofPixels& pix);
    const ofPixels& out() const;
    void setThreshold(float threshold);
    float getThreshold() const;
    void setBinary(bool mode);
    bool isBinary() const;
    void setInvert(bool mode);
    bool isInvert() const;
};





} // end namespace
