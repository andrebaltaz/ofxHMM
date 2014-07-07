#include "testApp.h"

//--------------------------------------------------------------
void testApp::setup(){	 
	counter = 0;
	vagRounded.loadFont("vag.ttf", 32);
	ofBackground(50,50,50);
    trained=false;
    
     
}

//--------------------------------------------------------------
void testApp::update(){
	counter = counter + 0.033f;
}

//--------------------------------------------------------------
void testApp::draw(){
	
	sprintf (timeString, "time: %i:%i:%i \nelapsed time %lli", ofGetHours(), ofGetMinutes(), ofGetSeconds(), ofGetElapsedTimeMillis());
	
	float w = vagRounded.stringWidth(eventString);
	float h = vagRounded.stringHeight(eventString);
	
	ofSetHexColor(0xffffff);
	vagRounded.drawString(eventString, 98,198);
	
	ofSetColor(255,122,220);
	vagRounded.drawString(eventString, 100,200);
	
	
	ofSetHexColor(0xffffff);
	vagRounded.drawString(timeString, 98,98);
	
	ofSetColor(255,122,220);
	vagRounded.drawString(timeString, 100,100);
	
}


//--------------------------------------------------------------
void testApp::keyPressed  (int key){ 
	sprintf(eventString, "keyPressed = (%i)", key);
    
    switch (key) {
            case ' ':
            point= new testHMM(vec_data);
            trained=true;
            break;
        case '1':
            point->Run(1, vec_data);
             vec_data.clear();
            break;
        case '2':
            point->Run(2, vec_data);
             vec_data.clear();
            break;
        case '3':
            point->Run(3, vec_data);
             vec_data.clear();
            break;
        case '4':
            point->Run(1, test_data);
           // vec_data.clear();
            break;
       
    }
    
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){ 
	sprintf(eventString, "keyReleased = (%i)", key);	
}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){
	sprintf(eventString, "mouseMoved = (%i,%i)", x, y);
}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){
	sprintf(eventString, "mouseDragged = (%i,%i - button %i)", x, y, button);
    mouse_coords.set(x,y);
    mousepos.push_back(mouse_coords);
}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
	sprintf(eventString, "mousePressed = (%i,%i - button %i)", x, y, button);

    

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){
	sprintf(eventString, "mouseReleased = (%i,%i - button %i)", x, y, button);
    int word=0;
        
    for (int i=1; i<mousepos.size(); i++) {
        float result;
        result=atan2((-1*(mousepos[i].y-mousepos[i-1].y)),(mousepos[i].x-mousepos[i-1].x))* 180 / PI;
        
        if (result<0) {
            word=(result+360)/30;
        }
        else
            word=result/30;
        
        if (word_data.size()>=100){
            word_data.erase(word_data.begin()+0);
        }
        
        word_data.push_back(word);
         //printf("vec %d=%f %f word=%d\n", i, mousepos[i].x-mousepos[i-1].x, mousepos[i].y-mousepos[i-1].y, word_data[i]);
    }
        while (word_data.size()<=100) {
            word_data.push_back(word);
        }

        vec_data.push_back(word_data);
    
    
    
    if (trained==false)
    {
        
        test_data.push_back(word_data);
        
//        for (int i=0; i<20; i++) {
//            vec_data.push_back(word_data);
//        }
    }

    
    for (int i=0; i<word_data.size(); i++) {
        printf("%d ", word_data[i]);
    }
    printf("\n");
    printf("% ld %ld\n", word_data.size(), vec_data.size());
    word_data.clear();
    mousepos.clear();
  
    
}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){
	sprintf(eventString, "resized = (%i,%i)", w, h);
}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){
	sprintf(eventString, "gotMessage %s ", msg.message.c_str());
}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 
	sprintf(eventString, "%i files dragged into the window at (%i, %i)", (int)dragInfo.files.size(), (int)dragInfo.position.x, (int)dragInfo.position.y);
}
