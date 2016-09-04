#include "ofxPixelUtility.h"

// begin namespace
namespace ofxPixelUtility {

/// color functions:

//-------------------------------------------------------------

vector<float> getMeanColor(const ofPixels& pixels){
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
    }

    const int numChannels = pixels.getNumChannels();
    const uint32_t length = pixels.getWidth()*pixels.getHeight();
    vector<float> result(numChannels);

    const unsigned char * pix = pixels.getData();

    for (int k = 0; k < numChannels; ++k){
        float sum = 0;
        for (uint32_t i = 0; i < length; ++i){
            sum += pix[i*numChannels+k];
        }
        result[k] = sum / (length*255.f);
    }

    return result;
}

//----------------------------------------------------------

float getMeanColor(const ofPixels& pixels, int channel){
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
    }

    const int numChannels = pixels.getNumChannels();
    channel = min(channel, numChannels-1);
    const uint32_t length = pixels.getWidth()*pixels.getHeight();
    
	float sum = 0.f;
    const unsigned char * pix = pixels.getData();
	
    for (uint32_t i = 0; i < length; ++i){
        sum += pix[i*numChannels+channel];
    }
    sum = sum / (length*255.f);

    return sum;
}

//----------------------------------------------------------

float getBrightness(const ofPixels& pixels){
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
    }

    const uint32_t length = pixels.getWidth()*pixels.getHeight()*pixels.getNumChannels();
    float sum = 0.f;
   
    const unsigned char* pix = pixels.getData();

    for (uint32_t i = 0; i < length; ++i){
        sum += pix[i];
    }

    sum = sum / (length*255.f);

    return sum;
}

//---------------------------------------------------------------

void adjustGain(ofPixels &pixels, float gain){
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
        return;
    }

    if (gain == 1.f){
        return;
    }

    gain = std::max(0.f, gain);

    const int32_t length = pixels.getWidth()*pixels.getHeight()*pixels.getNumChannels();
    unsigned char * pix = pixels.getData();

    for (uint32_t i = 0; i < length; ++i){
        float temp = (pix[i] / 255.f) * gain;
        temp = std::min(1.f, temp);
        pix[i] = static_cast<int>(temp * 255.f + 0.5f);
    }
}

void adjustGain(ofPixels &pixels, float gain, int channel){
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
        return;
    }

    if (gain == 1.f){
        return;
    }

    gain = std::max(0.f, gain);

    const int numChannels = pixels.getNumChannels();
    channel = std::max(0, std::min(numChannels-1, channel));
    const int32_t length = pixels.getWidth()*pixels.getHeight();

    unsigned char * pix = pixels.getData();

    for (int i = 0; i < length; ++i){
        int index = i * numChannels + channel;
        float temp = (pix[index] / 255.f) * gain;
        temp = std::min(1.f, temp);
        pix[index] = static_cast<int>(temp * 255.f + 0.5f);
    }
}

//------------------------------------------------------------------------------------

void adjustBrightnessAndContrast(ofPixels& pixels, float brightness, float contrast){
    if (pixels.isAllocated()){
        const uint32_t length = pixels.getWidth()*pixels.getHeight()*pixels.getNumChannels();

        contrast = max(0.f, contrast);

        unsigned char * pix = pixels.getData();

        // only contrast:
        if ((brightness == 0)&&(contrast != 1)){
            for (uint32_t i = 0; i < length; ++i){
                // simple formula for contrast adjustment
                int temp = static_cast<int>((pix[i] - 128)*contrast + 128.5f);
                temp = max(0, min(255, temp));
                pix[i] = static_cast<unsigned char>(temp);
            }
        }
        // only brightness:
        else if ((brightness != 0)&&(contrast == 1)){
            for (uint32_t i = 0; i < length; ++i){
                int temp = static_cast<int>(pix[i] + brightness);
                temp = max(0, min(255, temp));
                pix[i] = static_cast<unsigned char>(temp);
            }
        }
        // brightness + contrast:
        else if ((brightness != 0)&&(contrast != 1)){
            for (uint32_t i = 0; i < length; ++i){
                int temp = static_cast<int>((pix[i] - 128)*contrast + 128.5 + brightness);
                temp = max(0, min(255, temp));
                pix[i] = static_cast<unsigned char>(temp);
            }
        // neither brightness nor contrast - do nothing.
        }
    } else {
        cout << "pixels not allocated!\n";
    }
}

void adjustBrightnessAndContrast(ofPixels& pixels, float brightness, float contrast, int channel){
    if (pixels.isAllocated()){
        const int numChannels = pixels.getNumChannels();
        const uint32_t length = pixels.getWidth()*pixels.getHeight();

        channel = max(0, min(numChannels-1, channel));
        contrast = max(0.f, contrast);

        unsigned char * pix = pixels.getData();

        // only contrast:
        if ((brightness == 0)&&(contrast != 1)){
            for (uint32_t i = 0; i < length; ++i){
                int index = i * numChannels + channel;
                // only contrast adjustment
                int temp = static_cast<int>((pix[index] - 128)*contrast + 128.5f);
                temp = max(0, min(255, temp));
                pix[index] = static_cast<unsigned char>(temp);
            }
        }
        // only brightness:
        else if ((brightness != 0)&&(contrast == 1)){
            for (uint32_t i = 0; i < length; ++i){
                int index = i * numChannels + channel;
                // only brightness adjustment
                int temp = static_cast<int>(pix[index] + brightness);
                temp = max(0, min(255, temp));
                pix[index] = static_cast<unsigned char>(temp);
            }
        }
        // brightness + contrast:
        else if ((brightness != 0)&&(contrast != 1)){
            for (uint32_t i = 0; i < length; ++i){
                int index = i * numChannels + channel;
                // simple formula for contrast and brightness adjustment
                int temp = static_cast<int>((pix[i] - 128)*contrast + 128.5f + brightness);
                temp = max(0, min(255, temp));
                pix[index] = static_cast<unsigned char>(temp);
            }
        // neither brightness nor contrast - do nothing
        }
    } else {
        cout << "pixels not allocated!\n";
    }

}

//---------------------------------------------------------------------------

void adjustHSB(ofPixels& pixels, float hue, float saturation, float brightness){
    if (hue == 0 && saturation == 1 && brightness == 1){
        // no adjustment
        return;
    }
    if (!pixels.isAllocated()){
        return;
    }

    const int channels = pixels.getNumChannels();

    if (channels == 3 || channels == 4){
        const uint32_t length = pixels.getWidth()*pixels.getHeight()*channels;
        unsigned char * pix = pixels.getData();

        for (uint32_t i = 0; i < length; i += channels){
            int r = pix[i];
            int g = pix[i+1];
            int b = pix[i+2];
            float newHue, newSaturation, newBrightness;
            // first calculate orginal HSB values:
            // get strongest channels
            float maxVal = r;
            if(g > maxVal) {
                maxVal = g;
            }
            if(b > maxVal) {
                maxVal = b;
            }
            // get weakest channel
            float minVal = r;
            if(g < minVal) {
                minVal = g;
            }
            if(b < minVal) {
                minVal = b;
            }
            // check for grays
            if(maxVal == minVal) {
                newHue = 0.f;
                newSaturation = 0.f;
                newBrightness = maxVal;
            } else {
                float hueSixth;
                if(r == maxVal) {
                    hueSixth = (g - b) / (maxVal - minVal);
                    if(hueSixth < 0.f)
                        hueSixth += 6.f;
                } else if (g == maxVal) {
                    hueSixth = 2.f + (b - r) / (maxVal - minVal);
                } else {
                    hueSixth = 4.f + (r - g) / (maxVal - minVal);
                }
                newHue = 255.f * hueSixth / 6.f;
                newSaturation = 255.f * (maxVal - minVal) / maxVal;
                newBrightness = maxVal;
            }
            //
            // adjust HSB:
            //
            newHue += hue*255.f;
            newSaturation *= saturation;
            newBrightness *= brightness;
            // wrap/clip new values
            newHue = fmodf(newHue, 255.f);
            if (newHue < 0) {
                newHue += 255.f;
            }
            newSaturation = max(0.f, min(255.f, newSaturation));
            newBrightness = max(0.f, min(255.f, newBrightness));
            //
            // convert back to RGB
            //
            if(newBrightness == 0) { // black
                pix[i] = 0;
                pix[i+1] = 0;
                pix[i+2] = 0;
            } else if(newSaturation == 0) { // grays
                pix[i] = newBrightness;
                pix[i+1] = newBrightness;
                pix[i+2] = newBrightness;
            } else {
                float hueSix = newHue * 6.f / 255.f;
                float saturationNorm = newSaturation / 255.f;
                int hueSixCategory = (int) floorf(hueSix);
                float hueSixRemainder = hueSix - hueSixCategory;
                unsigned char pv = (unsigned char) ((1.f - saturationNorm) * newBrightness);
                unsigned char qv = (unsigned char) ((1.f - saturationNorm * hueSixRemainder) * newBrightness);
                unsigned char tv = (unsigned char) ((1.f - saturationNorm * (1.f - hueSixRemainder)) * newBrightness);
                switch(hueSixCategory) {
                    case 0: case 6: // r
                        pix[i] = newBrightness;
                        pix[i+1] = tv;
                        pix[i+2] = pv;
                        break;
                    case 1: // g
                        pix[i] = qv;
                        pix[i+1] = newBrightness;
                        pix[i+2] = pv;
                        break;
                    case 2:
                        pix[i] = pv;
                        pix[i+1] = newBrightness;
                        pix[i+2] = tv;
                        break;
                    case 3: // b
                        pix[i] = pv;
                        pix[i+1] = qv;
                        pix[i+2] = newBrightness;
                        break;
                    case 4:
                        pix[i] = tv;
                        pix[i+1] = pv;
                        pix[i+2] = newBrightness;
                        break;
                    case 5: // back to r
                        pix[i] = newBrightness;
                        pix[i+1] = pv;
                        pix[i+2] = qv;
                        break;
                }
            }
        }
    } else {
        cout << "Error: pixels must be RGB or RGBA\n";
        return;
    }
}


//-----------------------------------------------------------------------------

void applyColorCurve(ofPixels& pixels, vector<unsigned char>& table, int channel){
    if (table.size() < 256){
        cout << "look up table size must not be smaller than 256!\n";
        return;
    }
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
        return;
    }

    int numChannels = pixels.getNumChannels();
    const uint32_t length = pixels.getWidth()*pixels.getHeight();
    channel = max(0, min(numChannels-1, channel));

    unsigned char * pix = pixels.getData();

    for (uint32_t i = 0; i < length; ++i){
        int index = i * numChannels + channel;
        // get color
        int color = pix[index];
        // look up color in table and assign new value
        pix[index] = table[color];
    }
}

void applyColorCurve(ofPixels& pixels, vector<unsigned char>& table){
    if (table.size() < 256){
        cout << "look up table size must not be smaller than 256!\n";
        return;
    }
    if (!pixels.isAllocated()){
        cout << "pixels not allocated!\n";
        return;
    }

    uint32_t length = pixels.getWidth()*pixels.getHeight()*pixels.getNumChannels();

    unsigned char * pix = pixels.getData();

    for (uint32_t i = 0; i < length; ++i){
        // get color
        int color = pix[i];
        // look up color in table and assign new value
        pix[i] = table[i];
    }
}

//------------------------------------------------------------------------

void convertToGrayscale(const ofPixels& src, ofPixels& dest, GrayscaleMode mode){
    if (!src.isAllocated()){
        cout << "source pixels not allocated!\n";
        return;
    }

    const int w = src.getWidth();
    const int h = src.getHeight();
    const int channels = src.getNumChannels();
    const uint32_t length = w*h*channels;
    
	if (channels != 3 && channels != 4){
        cout << "source pixels must have 3 or 4 channels!\n";
		return;
    }

    if (w != dest.getWidth() || h != dest.getHeight() || dest.getNumChannels() != 1){
        dest.allocate(w, h, 1);
    }

    const unsigned char * pix = src.getData();
    unsigned char * out = dest.getData();

	switch (mode){
        case GrayscaleMode::LIGHTNESS:
            for (int i = 0, k = 0; i < length; i+=channels, k++){
                float sum = 0;
                // The formula for luminosity is 0.21 R + 0.72 G + 0.07 B.
                int r = pix[i];
                int g = pix[i+1];
                int b = pix[i+2];

                unsigned int min = r;
                if (min > g) min = g;
                if (min > b) min = b;

                unsigned int max = r;
                if (max > g) max = g;
                if (max > b) max = b;

                out[k] = (float)(min+max)/2.f;
			}
		break;
        case GrayscaleMode::AVERAGE:
            for (int i = 0, k = 0; i < length; i+=channels, k++){
                float sum = 0;
                sum += pix[i];
                sum += pix[i+1];
                sum += pix[i+2];
                out[k] = sum/3.f;
			}
		break;
        case GrayscaleMode::LUMINANCE:
            for (int i = 0, k = 0; i < length; i+=channels, k++){
                float sum = 0;
                // The formula for luminosity is 0.21 R + 0.72 G + 0.07 B.
                sum += pix[i]*0.21f;
                sum += pix[i+1]*0.72f;
                sum += pix[i+2]*0.07f;
                out[k] = sum;
			}
		break;
	}
}

//-------------------------------------------------------

void applyAlphaMask(const ofPixels& src, const ofPixels& mask, ofPixels& dest, bool invertMask){
    const int channels = src.getNumChannels();

    if (channels != 1 && channels != 3){
		cout << "source image must be RGB or greyscale!\n";
        return;
	}

    if (mask.getNumChannels() != 1){
		cout << "mask must be greyscale!\n";
        return;
	}
	
    const int w = src.getWidth();
    const int h = src.getHeight();

    if (w != mask.getWidth() || h != mask.getHeight()){
		"image and mask dimensions don't match!\n";
        return;
	}
	
    if (dest.getWidth() != w || dest.getHeight() != h || dest.getNumChannels() != (channels+1)){
        dest.allocate(w, h, channels+1);
    }

    const unsigned char* in = src.getData();
    const unsigned char* m = mask.getData();
    unsigned char* out = dest.getData();
    const uint32_t size = w*h;

    // source image and mask are both 1 channel (G), result is 2 channel (GA)
    if (channels == 1){
        for (uint32_t i = 0; i < size; ++i){
            const int k = 2*i;
            out[k] = in[i];
            out[k+1] = invertMask ? 255 - m[i] : m[i];
        }
    }
    // source image is 3 channel (RGB) and mask is 1 channel (G), result is 4 channel (RGBA)
    else {
        for (uint32_t i = 0; i < size; ++i){
            const int j = 3*i;
            const int k = 4*i;
            out[k] = in[j];
            out[k+1] = in[j+1];
            out[k+2] = in[j+2];
            out[k+3] = invertMask ? 255 - m[i] : m[i];
        }
    }
}

//-------------------------------------------------------

/// Movement Detection

MovementDetection::MovementDetection() {
    /// default constructor
    myFeedback = 0.f;
    myGain = 1.f;
    myThreshold = 0;
    bBinary = false;
    bInvert = false;
    width = 0;
    height = 0;
    numChannels = 0;
    bAllocated = false;
    buf_1 = nullptr;
    buf_w = nullptr;
}

MovementDetection::MovementDetection(int w, int h, int channels) {
    myFeedback = 0.f;
    myGain = 1.f;
    myThreshold = 0;
    bBinary = false;
    bInvert = false;
    width = 0;
    height = 0;
    numChannels = 0;
    bAllocated = false;
    buf_1 = nullptr;
    buf_w = nullptr;

    allocate(w, h, channels);
}

MovementDetection::~MovementDetection() {
    /// destructor which deletes the buffers
    if (bAllocated) {
        delete[] buf_1;
        delete[] buf_w;
    }
}

void MovementDetection::allocate(int w, int h, int channels){
    w = std::max(0, w);
    h = std::max(0, h);
    channels = std::max(0, channels);

    uint32_t size = w * h * channels;

    // check for bad dimensions
    if (size == 0 || channels > 4){
        cout << "bad dimensions, not allocating!\n";
        return;
    }

    if (bAllocated){
        delete[] buf_1;
        delete[] buf_w;
    }

    buf_1 = new float[size];
    buf_w = new float[size];

    outPixels.allocate(w, h, channels);

    width = w;
    height = h;
    numChannels = channels;
    bAllocated = true;

    clearFilter();
}


bool MovementDetection::isAllocated(){
    return bAllocated;
}

void MovementDetection::in(const ofPixels& inPixels) {
    if (!inPixels.isAllocated()){
        cout << "incoming pixels not allocated!\n";
        return;
    }

    const int w = inPixels.getWidth();
    const int h = inPixels.getHeight();
    const int channels = inPixels.getNumChannels();

    /// check if the dimensions have changed and reallocate if necessary.

    if ((w != width)||(h != height)||(channels != numChannels)){
        allocate(w, h, channels);
    }

    if (!bAllocated){
        cout << "error: allocation failed!\n"; // shouldn't happen!
        return;
    }

    /// calculate the following difference equation:
    /// 
	/// norm = (1-fb)*gain
    /// w[n] = x[n]*norm + w[n-1]*fb
    /// y[n] = abs(w[n] - w[n-1])


    const unsigned char* inPix = inPixels.getData();
    unsigned char* outPix = outPixels.getData();
    const uint32_t size = w*h*channels;
    const float norm = (1.f - myFeedback);
    const float gain = myGain * 255.f;
	// generate tiny random offset to protect against denormals
    const float noise = ofRandom(-1e-006, 1e-006);

    // special case 1: everything is low
    if (myThreshold > 255){
        const int low = (bBinary && bInvert) ? 255 : 0;

        for (uint32_t i = 0; i < size; ++i){
            buf_w[i] = inPix[i] / 255.f; // update the buffer nevertheless
            outPix[i] = low;
        }
    }
    // special case 2: everything is high
    else if (bBinary && myThreshold == 0){
        const int high = bInvert ? 0 : 255;

        for (uint32_t i = 0; i < size; ++i){
            buf_w[i] = inPix[i] / 255.f; // update the buffer nevertheless
            outPix[i] = high;
        }
    }
    /*
    // special case 3: no need for threshold checking
    else if (!bBinary && myThreshold == 0){
        for (uint32_t i = 0; i < size; ++i){
            // averaging
            buf_w[i] = inPix[i] * norm / 255.f + buf_1[i] * myFeedback + noise;
            // take difference, amplify and round to int
            int32_t temp = static_cast<int32_t>((buf_w[i] - buf_1[i])*gain + 0.5f);
            // take absolute value
            temp = abs(temp);
            // clip high output
            temp = std::min(255, temp);
            outPix[i] = static_cast<unsigned char>(temp);
        }
    }
    */
    // binary thresholding
    else if (bBinary){
        const int low = bInvert ? 255 : 0;
        const int high = bInvert ? 0 : 255;

        for (uint32_t i = 0; i < size; ++i){
            // averaging
            buf_w[i] = inPix[i] * norm / 255.f + buf_1[i] * myFeedback + noise;
            // take difference, amplify and round to int
            int32_t temp = static_cast<int32_t>((buf_w[i] - buf_1[i])*gain + 0.5f);
            // take absolute value
            temp = abs(temp);
            // apply threshold function
            temp = (temp >= myThreshold) ? high : low;
            outPix[i] = static_cast<unsigned char>(temp);
        }
    }
    // non-binary thresholding
    else {
        for (uint32_t i = 0; i < size; ++i){
            // averaging
            buf_w[i] = inPix[i] * norm / 255.f + buf_1[i] * myFeedback + noise;
            // take difference, amplify and round to int
            int32_t temp = static_cast<int32_t>((buf_w[i] - buf_1[i])*gain + 0.5f);
            // take absolute value
            temp = abs(temp);
            // apply threshold function and clip high output
            if (temp >= myThreshold){
                temp = std::min(255, temp);
            } else {
                temp = 0;
            }
            outPix[i] = static_cast<unsigned char>(temp);
        }
    }

    /// move buffers by swapping the pointers.
    float * temp = buf_1;
    buf_1 = buf_w;
    buf_w = temp;
    // buf_w points now to old buf_1 which will be overwritten the next time.
}

const ofPixels& MovementDetection::out() const {
    /// returns ofPixels (don't need to check. in the worst case, outPixels is simply not allocated.
    return outPixels;
}

void MovementDetection::clearFilter() {
    /// clear the buffers (but only if the filter is not empty)
    const uint32_t size = width*height*numChannels;
    if (bAllocated){
        for (uint32_t i = 0; i < size; ++i){
            buf_1[i] = 0;
        }
    }
}

void MovementDetection::setAveraging(float coeff) {
    myFeedback = max(0.f, min(0.999999f, coeff));
}

void MovementDetection::setGain(float gain) {
    myGain = gain;
}

void MovementDetection::setBinary(bool mode) {
    bBinary = mode;
}

void MovementDetection::setThreshold(float threshold) {
    // 0.f = everything passes, 1.f = nothing passes
    myThreshold = max(0, static_cast<int>(threshold*256.f + 0.5f));
}

void MovementDetection::setInvert(bool mode) {
    bInvert = mode;
}

float MovementDetection::getAveraging() const {
    return myFeedback;
}

float MovementDetection::getGain() const {
    return myGain;
}

bool MovementDetection::isBinary() const {
    return bBinary;
}

float MovementDetection::getThreshold() const {
    return (float) myThreshold / 256.f;
}

bool MovementDetection::isInvert() const {
    return bInvert;
}


//---------------------------------------------------

BackgroundSubstraction::BackgroundSubstraction(){
    myThreshold = 0;
    bBinary = false;
    bInvert = false;
}

BackgroundSubstraction::~BackgroundSubstraction(){

}

void BackgroundSubstraction::setBackground(const ofPixels &pix){
    if (pix.isAllocated()){
        myBackground = pix;
        width = pix.getWidth();
        height = pix.getHeight();
        numChannels = pix.getNumChannels();
        outPixels.allocate(width, height, numChannels);
    }
}

void BackgroundSubstraction::setBackground(int w, int h, int channels, ofColor color){
    myBackground.allocate(w, h, channels);
    myBackground.setColor(color);
    width = w;
    height = h;
    numChannels = channels;
    outPixels.allocate(width, height, numChannels);
}

const ofPixels& BackgroundSubstraction::getBackground() const {
    return myBackground;
}

bool BackgroundSubstraction::isBackgroundSet() const {
    return myBackground.isAllocated();
}

void BackgroundSubstraction::in(const ofPixels &pix){
    if (pix.getWidth() != width || pix.getHeight() != height || pix.getNumChannels() != numChannels){
        cout << "input doesn't match background dimensions!\n";
        return;
    }

    if (!isBackgroundSet()){
        cout << "background is not set!\n";
        return;
    }

    const unsigned char* in = pix.getData();
    const unsigned char* back = myBackground.getData();
    unsigned char* out = outPixels.getData();
    const uint32_t size = width * height * numChannels;

    if (!bBinary){
        for (uint32_t i = 0; i < size; ++i){
            int temp = (int) in[i] - (int) back[i];
            temp = abs(temp);
            out[i] = (temp >= myThreshold) ? temp : 0;
        }
    } else {
        const int low = bInvert ? 255 : 0;
        const int high = bInvert ? 0 : 255;

        for (uint32_t i = 0; i < size; ++i){
            int temp = (int) in[i] - (int) back[i];
            temp = abs(temp);
            temp = (temp >= myThreshold) ? high : low;
            out[i] = temp;
        }
    }
}


const ofPixels& BackgroundSubstraction::out() const{
    return outPixels;
}

void BackgroundSubstraction::setThreshold(float threshold){
    // 0.f = everything passes, 1.f = nothing passes
    myThreshold = std::max(0, static_cast<int>(threshold * 256.f + 0.5f));
}

float BackgroundSubstraction::getThreshold() const {
    return myThreshold / 256.f;
}


void BackgroundSubstraction::setBinary(bool mode){
    bBinary = mode;
}

bool BackgroundSubstraction::isBinary() const {
    return bBinary;
}

void BackgroundSubstraction::setInvert(bool mode){
    bInvert = mode;
}

bool BackgroundSubstraction::isInvert() const {
   return bInvert;
}


} // end namespace
