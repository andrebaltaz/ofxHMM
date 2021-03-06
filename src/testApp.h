#pragma once

#include "ofMain.h"
#include "HMM.h"

class testApp : public ofBaseApp{
	
public:
    
    void setup();
    void update();
    void draw();
    
    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    float 			counter;
    ofTrueTypeFont 	vagRounded;
    char eventString[255];
    char timeString[255];
    ofVec2f mouse_coords;
    vector<ofVec2f> mousepos;
    
    vector<int> word_data;
    vector<vector<int> >vec_data;
    vector<vector<int> >test_data;
    bool trained;
    testHMM* point;
};

