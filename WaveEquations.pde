import pallav.Matrix.*;

import peasy.*; //from a libary called PeasyCam so that we can move the camera
PeasyCam p;
float bohr = 1.5;
float bondl = 1.5;
PVector [][][] SI;
ArrayList<PVector> [] SI3;
int tz = 1;
PVector c1 = new PVector(0, 0, -bondl);
PVector c2 = new PVector(0, 0, bondl);
PVector h1 = new PVector(0,1.73,2.5);
PVector h2 = new PVector(0,1.73,-2.5);
PVector h3 = new PVector(0,-1.73,2.5);
PVector h4 = new PVector(0,-1.73,-2.5);
PVector h;
float dlb = -6;
float dub = 6;
float step = .1;
float [] pSet;
int loops= 0;
int shells;
final int inputStep = 100;
final int inputStep2 = 4;
void setup() {
  //noFill();
  SI3 = new ArrayList[inputStep];
  size(800, 800, P3D);
  p = new PeasyCam(this, 400, 400, 400, 400);
  colorMode(HSB);
  SI = new PVector[width/inputStep2][width/inputStep][15];
}

void draw() {
  pSet = new float[] {0.003, 0.0110, 0.005, 0.005, 0.005, 0.005, 0.0100, 0.0100, 0.0100, 0.0100, 0.00001, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 0.0100, 
    0.0100, 0.0100, 0.0100, 0.0100};
  shells = 4;
  background(0);
  drawAxes();
  translate(400, 400, 400);
  //sketchParam(new PVector(0,0,0),new PVector(0,0),-2.0,2.0,0,TWO_PI,.1,PI/5);
  //sketchDensityWave();
  if (loops == 0) {
    //makeDensityGraphRadialA();
    makeDensityGraph3D();
    loops++;
  }
  //sketchParamGraph();
  sketchParamGraph3D();
  ArrayList<PVector> h = new ArrayList<PVector>();
  h.add(c1);
  h.add(c2);
  h.add(h1);
  h.add(h2);
  h.add(h3);
  h.add(h4);
  strokeWeight(20);

  for(int i = 0; i < h.size(); i++) {
    stroke(120,0,255);
    if (i <= 1) {
      stroke(80,255,255);}
    point(m2(h.get(i).x),m2(h.get(i).y),m2(h.get(i).z));
  }
}

void makeDensityGraphRadialA() {
  strokeWeight(11);
  int z1 = 6;
  int z2 = 6;
  for (int shell = 3; shell < 1 + shells; shell++) {
    float pSet1 = pSet[shell];
    for (int z = 0; z < SI.length; z++) {
      float rad = dub;
      while (lcao12(new PVector(rad, 0, map(z, 0, SI.length, dlb, dub)),6,6,4) - pSet1 < 0) {
        rad-= .001;
        if (rad < 0) {
          rad = 0; 
          break;}}
      for (int theta = 0; theta < SI[0].length; theta++) {
        float thetaReal = map(theta, 0, SI[0].length, 0, TWO_PI);
        SI[z][theta][shell-1] = new PVector(rad*cos(thetaReal), rad*sin(thetaReal), map(z, 0, SI.length, dlb, dub));
      }}}
  strokeWeight(2);
}

void makeDensityGraph3D() {
  strokeWeight(2);
  float pThresh = .005;
float step = .1;
  for (int x = 0; x < SI3.length; x++) {
    ArrayList<PVector> points = new ArrayList<PVector>();
    int old = -1;
    for (float y = dlb; y < dub; y+=step) {
      for (float z = dlb; z < dub; z+=step) {
        float t = lcao12(new PVector(map(x,0,SI3.length,dlb,dub),y,z), 6,6, 4)-pThresh;
        int thresh = int(t/abs(t));
        if (thresh != old) {
          points.add(new PVector(map(x,0,SI3.length,dlb,dub), y, z));}
          old = thresh;}}
    //int old = -1;
    //for (float y = dlb; y < 0; y+=ms) {
    //  old = -1;
    //  for (float x = dlb; x < 0; x+=ms) {
    //    float t = hybridBond(new PVector(x,y,map(z,0,SI3.length,dlb,dub)), 6,6, 4)-pThresh;
    //    int thresh = int(t/abs(t));
    //    if (thresh != old) {
    //      points.add(new PVector(x,y,map(z,0,SI3.length,dlb,dub)));
    //      old = -1;
    //      break;}}}
    //old = -1;
    //for (float y = dub; y > 0; y-=ms) {
    //  old = -1;
    //  for (float x = dlb; x < 0; x+=ms) {
    //    float t = hybridBond(new PVector(x,y,map(z,0,SI3.length,dlb,dub)), 6,6, 4)-pThresh;
    //    int thresh = int(t/abs(t));
    //    if (thresh != old) {
    //      points.add(new PVector(x,y,map(z,0,SI3.length,dlb,dub)));
    //      old = -1;
    //      break;}}}
    //old = -1;
    //for (float y = dlb; y < 0; y+=ms) {
    //  old = -1;
    //  for (float x = dub ; x >0; x-=ms) {
    //    float t = hybridBond(new PVector(x,y,map(z,0,SI3.length,dlb,dub)), 6,6, 4)-pThresh;
    //    int thresh = int(t/abs(t));
    //    if (thresh != old) {
    //      points.add(new PVector(x,y,map(z,0,SI3.length,dlb,dub)));
    //      old = -1;
    //      break;}}}
    //old = -1;
    //for (float y = dub; y > 0; y-=ms) {
    //  old = -1;
    //  for (float x = dub ; x >0; x-=ms) {
    //    float t = hybridBond(new PVector(x,y,map(z,0,SI3.length,dlb,dub)), 6,6, 4)-pThresh;
    //    int thresh = int(t/abs(t));
    //    if (thresh != old) {
    //      points.add(new PVector(x,y,map(z,0,SI3.length,dlb,dub)));
    //      old = -1;
    //      break;}}}
    SI3[x] = distanceSort(points);
  }
}
ArrayList<PVector> distanceSort(ArrayList<PVector> in) {
    ArrayList<PVector> n = new ArrayList<PVector>();
    boolean [] done = new boolean[in.size()];
    int index = 0;
    for (int i = 0; i < in.size(); i++) {
      if (in.get(i).x < in.get(index).x) {
        index = i;}}
    if (in.size() == 0) {return n;}
    PVector o = in.get(index);
    done[index] = true;
    n.add(new PVector(o.x,o.y,o.z));
    for (int i = 0; i < in.size(); i++) {
      float minDist = 100;
      int i1 = 0;
      for (int j = 0; j < in.size(); j++) {
        if (i != 0) {
          if (!done[j]) {
            if (PVector.dist(n.get(i),in.get(j))<minDist) {
              i1 = j;
              minDist = PVector.dist(n.get(i),in.get(j));}}}
        else {
          if (!done[j] && in.get(j).y > n.get(0).y) {
            if (PVector.dist(n.get(i),in.get(j))<minDist) {
              i1 = j;
              minDist = PVector.dist(n.get(i),in.get(j));}}}      
        }
      done[i1] = true;
      o = in.get(i1);
      n.add(new PVector(o.x,o.y,o.z));}
    return n;}

void sketchParamGraph3D() {
  strokeWeight(2);
  for (int i = 1; i < SI3.length; i++) {
    for (int j = 0; j < SI3[i].size(); j++) {
      int jMod = (j+1+SI3[i].size())%SI3[i].size();
      line(m2(SI3[i].get(j).x),m2(SI3[i].get(j).y),m2(SI3[i].get(j).z),
           m2(SI3[i].get(jMod).x),m2(SI3[i].get(jMod).y),m2(SI3[i].get(jMod).z));
      float chunk = float(SI3[i].size())/SI3[i-1].size();
      for (int k = 0; k < SI3[i-1].size(); k++) {
        jMod = int(chunk*k);
        line(m2(SI3[i].get(jMod).x),m2(SI3[i].get(jMod).y),m2(SI3[i].get(jMod).z),
             m2(SI3[i-1].get(k).x),m2(SI3[i-1].get(k).y),m2(SI3[i-1].get(k).z));
      }
    }
  }
}

void sketchParamGraph() {
  for (int shell =3; shell < shells; shell++) {
    for (int z = 1; z < SI.length; z++) {
      if (SI[z][0][shell-1].x != 0  && SI[z-1][0][shell-1].x != 0) {
        for (int theta = 0; theta < SI[0].length; theta++) {
          stroke(85, 255, 255);
          strokeWeight(2);
          int thetaModUp = (theta +1) % SI[0].length;
          //int thetaModDown = (theta + SI[0].length -1) % SI[0].length;
          line(m2(SI[z][theta][shell-1].x), m2(SI[z][theta][shell-1].y), m2(SI[z][theta][shell-1].z), 
            m2(SI[z-1][theta][shell-1].x), m2(SI[z-1][theta][shell-1].y), m2(SI[z-1][theta][shell-1].z));
          line(m2(SI[z][theta][shell-1].x), m2(SI[z][theta][shell-1].y), m2(SI[z][theta][shell-1].z), 
            m2(SI[z-1][thetaModUp][shell-1].x), m2(SI[z-1][thetaModUp][shell-1].y), m2(SI[z-1][thetaModUp][shell-1].z));
          //line(m2(SI[z][theta][shell-1].x), m2(SI[z][theta][shell-1].y), m2(SI[z][theta][shell-1].z), 
          //  m2(SI[z-1][thetaModDown][shell-1].x), m2(SI[z-1][thetaModDown][shell-1].y), m2(SI[z-1][thetaModDown][shell-1].z));
        }
      }
    }
  }
}

void sketchDensityWave() {
  strokeWeight(11);
  int z1 = 3;
  int z2 = 3;
  int o = 3;
  int power = 5;
  int num = int(pow(10, power));
  float thresh = 0.01;
  float delta = 0.001;
  for (int i = 0; i < num; i++) {
    PVector r = new PVector(dlb+random(dub-dlb), dlb+random(dub-dlb), dlb+random(dub-dlb));
    float p = lcao12(r, z1, z2, o);
    //getNormal(r, z1, z2, o);
    //println(p);
    if (p < thresh+delta && p > thresh-delta) {
      stroke(map(log(p), -12, -1, 0, 255), 255, 255);
      point(m2(r.x), m2(r.y), m2(r.z));
    }
  }
  strokeWeight(2);
}

PVector getNormal(PVector r, int z1, int z2, int o) {
  float mx = (lcao12(new PVector(r.x+step, r.y, r.z), z1, z2, o)-lcao12(new PVector(r.x-step, r.y, r.z), z1, z2, o))/2/step;
  float my = (lcao12(new PVector(r.x, r.y+step, r.z), z1, z2, o)-lcao12(new PVector(r.x, r.y-step, r.z), z1, z2, o))/2/step;
  float mz = (lcao12(new PVector(r.x, r.y, r.z+step), z1, z2, o)-lcao12(new PVector(r.x, r.y, r.z-step), z1, z2, o))/2/step;
  float p = sqrt(sq(mx)+sq(my)+sq(mz));
  PVector n = new PVector(mx/p, my/p, mz/p);
  return n;
}

PVector intersectPlanes(PVector r1, PVector r2, PVector r3, PVector n1, PVector n2, PVector n3) {
  PVector l1 = crossP(n1, n2);
  PVector l2 = crossP(n1, n3);
  PVector lb1 = new PVector(0, 0, 0);
  float [] c = new float[] {PVector.dot(n1, r1), PVector.dot(n2, r2)};
  float[] yz = multiply(inverse(new float[][] {{n1.y, n1.z}, {n2.y, n2.z}}), c);
  lb1.y =yz[0];
  lb1.z =yz[1];
  return null;}

PVector crossP(PVector v1, PVector v2) {
  return new PVector(v1.y*v2.x-v1.x*v2.y, v1.z*v2.x-v1.x*v2.z, v1.x*v2.y-v1.y*v2.x);
}

float[][] inverse(float [][] in) {
  float d = 1/(in[1][1]*in[0][0]-in[1][0]*in[0][1]);
  return new float[][] {{d*in[1][1], d*-1*in[0][1]}, {d*-1*in[0][1], d*in[0][0]}};
}
float[] multiply(float[][] one, float[] two) {
  return new float[] {one[0][0]*two[0]+one[0][1]*two[1], 
    one[1][0]*two[0]+one[1][1]+two[1]};
}



float lcao12(PVector in, int z1, int z2, int o) {
  float d1 = PVector.sub(in, c1).mag();
  float d2 = PVector.sub(in, c2).mag();
  float theta = atan2(in.y-c1.y, in.x-c1.x);
  float phi1 = atan2(sqrt(sq(in.x-c1.x)+sq(in.y-c1.y)), in.z-c1.z);
  float phi2 = atan2(sqrt(sq(in.x-c2.x)+sq(in.y-c2.y)), in.z-c2.z);
  return waveTransform(o, o, 1, d1, d2, theta, phi1, phi2, z1, z2);
}
float SP21(PVector in, int z1) {
  float d1 = PVector.sub(in, c1).mag();
  float theta = atan2(in.y-c1.y, in.x-c1.x);
  float phi1 = atan2(sqrt(sq(in.x-c1.x)+sq(in.y-c1.y)), in.z-c1.z);
  float pSum = 0.00000000;
  pSum += (sqrt(.33)*wave(2, d1, z1, theta, phi1)+sqrt(.67)*wave(4,d1,z1,theta,phi1));
  pSum += (sqrt(.33)*wave(2, d1, z1, theta, phi1)-sqrt(.16)*wave(4,d1,z1,theta,phi1)+sqrt(.5)*wave(5,d1,z1,theta,phi1));
  pSum += (sqrt(.33)*wave(2, d1, z1, theta, phi1)-sqrt(.16)*wave(4,d1,z1,theta,phi1)-sqrt(.5)*wave(5,d1,z1,theta,phi1));
  return pSum;
}
float SP22(PVector in, int z1) {
  float d1 = PVector.sub(in, c2).mag();
  float theta = atan2(in.y-c2.y, in.x-c2.x);
  float phi1 = atan2(sqrt(sq(in.x-c2.x)+sq(in.y-c2.y)), in.z-c2.z);
  float pSum = 0.00000000;
  pSum += (sqrt(.33)*wave(2, d1, z1, theta, phi1)+sqrt(.67)*wave(4,d1,z1,theta,phi1));
  pSum += (sqrt(.33)*wave(2, d1, z1, theta, phi1)-sqrt(.16)*wave(4,d1,z1,theta,phi1)+sqrt(.5)*wave(5,d1,z1,theta,phi1));
  pSum += (sqrt(.33)*wave(2, d1, z1, theta, phi1)-sqrt(.16)*wave(4,d1,z1,theta,phi1)-sqrt(.5)*wave(5,d1,z1,theta,phi1));
  return pSum;
}

float hybridBond(PVector in, int z1, int z2, int o) {
  ArrayList<PVector> h = new ArrayList<PVector>();
  h.add(h1);
  h.add(h2);
  h.add(h3);
  h.add(h4);
  
  float pSum = 0;
  for(int i = 0; i < h.size(); i++) {
    float theta = atan2(in.y-h.get(i).y, in.x-h.get(i).x);
    float phi2 = atan2(sqrt(sq(in.x-h.get(i).x)+sq(in.y-h.get(i).y)), in.z-h.get(i).z);
    float d2 = PVector.sub(in, h.get(i)).mag();
    pSum+= wave(1, d2, z2, theta, phi2);
  }
  pSum+= SP22(in,z1)+SP21(in,z2);// (like in ethane)
  return sq(pSum);
}
  

float waveTransform(int s, int e, int p, float r1, float r2, float theta, float phi1, float phi2, int z1, int z2) {
  float pSum = 0;
  for (int i = s; i <= e; i++) {
    pSum+=wave(i, r1, z1, theta, phi1)+p*(wave(i, r2, z2, theta, phi2));
  }
  return sq(pSum);
}

float wave(int i, float x, int Z, float theta, float phi) {
  if (i == 1) {
    return oneS(x, Z);
  } else if (i == 2) {
    return twoS(x, Z);
  } else if (i == 3) {
    return twoPx(x, Z, theta, phi);
  } else if (i == 4) {
    return twoPy(x, Z, theta, phi);
  } else if (i == 5) {
    return twoPz(x, Z, theta, phi);
  } else if (i == 6) {
    return threeS(x, Z);
  } else if (i == 7) {
    return threePx(x, Z, theta, phi);
  } else if (i == 8) {
    return threePy(x, Z, theta, phi);
  } else if (i == 9) {
    return threePz(x, Z, theta, phi);
  } else if (i == 10) {
    return threeDxx(x, Z, theta, phi);
  } else if (i == 11) {
    return threeDxy(x, Z, theta, phi);
  } else if (i == 12) {
    return threeDyx(x, Z, theta, phi);
  } else if (i == 13) {
    return threeDzy(x, Z, theta, phi);
  } else if (i == 14) {
    return threeDyz(x, Z, theta, phi);
  } else {
    return 0;
  }
}


float oneS(float r, int Z) {
  float p = Z*r/bohr;
  return (1/sqrt(PI))*pow(Z/bohr, 1.5)*exp(-abs(p));}
float twoS(float r, int Z) {
  float p = Z*r/bohr;
  return (1/sqrt(32*PI))*pow(Z/bohr, 1.5)*(2-p)*exp(-abs(p/2));}
float twoPx(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/sqrt(32*PI)*pow(Z/bohr, 1.5)*p*exp(-p/2)*cos(phi);}
float twoPy(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/sqrt(64*PI)*pow(Z/bohr, 1.5)*p*exp(-p/2)*sin(phi)*cos(theta);}
float twoPz(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/sqrt(64*PI)*pow(Z/bohr, 1.5)*p*exp(-p/2)*sin(phi)*sin(theta);}
float threeS(float r, int Z) {
  float p = Z*r/bohr;
  return 1/81/sqrt(3*PI)*pow(Z/bohr, 1.5)*(27-18*p+2*sq(p))*exp(-p/3);
}
float threePx(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/81*sqrt(2/PI)*pow(Z/bohr, 1.5)*(6*r-sq(p))*exp(-p/3)*cos(phi);
}
float threePy(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/81/sqrt(PI)*pow(Z/bohr, 1.5)*(6*r-sq(p))*exp(-r/3)*sin(phi)*cos(theta);
}

float threePz(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/81/sqrt(PI)*pow(Z/bohr, 1.5)*(6*r-sq(p))*exp(-r/3)*sin(phi)*sin(theta);
}
float threeDxx(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/81/sqrt(6*PI)*pow(Z/bohr, 1.5)*sq(p)*exp(-p/3)*(3*sq(cos(phi))-1);
}

float threeDxy(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return float(1)/81/sqrt(PI)*pow(Z/bohr, 1.5)*sq(p)*exp(-p/3)*sin(phi)*cos(phi)*cos(theta);
}
float threeDyx(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/81/sqrt(PI)*pow(Z/bohr, 1.5)*sq(p)*exp(-p/3)*sin(phi)*cos(phi)*sin(theta);
}
float threeDzy(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/162/sqrt(PI)*pow(Z/bohr, 1.5)*sq(p)*exp(-p/3)*sq(sin(phi))*cos(theta*2);
}
float threeDyz(float r, int Z, float theta, float phi) {
  float p = Z*r/bohr;
  return 1/162/sqrt(PI)*pow(Z/bohr, 1.5)*sq(p)*exp(-p/3)*sq(sin(phi))*sin(2*theta);
}

float x(float p, float theta, float phi) {
  return p*cos(theta)*sin(phi);
}
float y(float p, float theta, float phi) {
  return p*sin(theta)*sin(phi);
}
float z(float p, float theta, float phi) {
  return p*cos(phi);
}
float p(float x, float y, float z) {
  return sqrt(sq(x)+sq(y)+sq(z));
}
float theta(float x, float y, float z) {
  return atan2(y, x);
}
float phi(float x, float y, float z) {
  return atan2(sqrt(x)+sqrt(y), z);
}

PVector r(float x, float theta) {
  return new PVector(x, sq(wave(1, x, tz, 0, 0))*sin(theta), sq(wave(1, x, tz, 0, 0))*cos(theta));
}



void sketchParam(PVector orient, PVector rotate, float lb1, float ub1, float lb2, float ub2, float step1, float step2) {
  translate(orient.x, orient.y, orient.z);
  rotateZ(rotate.x);
  rotateX(rotate.y);
  for (float u = lb1; u < ub1; u+=step1) {
    for (float v = lb2; v < ub2; v+= step2) {
      if (u-step1>lb1) {
        stroke(255, 0, 255);
        strokeWeight(2);
        line(repack(r(u, v)).x, repack(r(u, v)).y, repack(r(u, v)).z
          , repack(r(u-step1, v)).x, repack(r(u-step1, v)).y, repack(r(u-step1, v)).z);
        line(repack(r(u, v)).x, repack(r(u, v)).y, repack(r(u, v)).z
          , repack(r(u-step1, v-step2)).x, repack(r(u-step1, v-step2)).y, repack(r(u-step1, v-step2)).z);
        line(repack(r(u, v)).x, repack(r(u, v)).y, repack(r(u, v)).z
          , repack(r(u-step1, v+step2)).x, repack(r(u-step1, v+step2)).y, repack(r(u-step1, v+step2)).z);
      }
    }
  }
}
PVector repack(PVector in) {
  return new  PVector(m1(in.x), m1(in.y), m1(in.z));
}

float m1(float in) {
  return map(in, 0, 4, 0, 400);
}

float m2(float in) {
  return map(in, dlb, dub, -400, 400);
}

void drawAxes() {
  strokeWeight(2);
  pushMatrix();
  translate(0, 0, 400);
  stroke(160, 255, 255);
  line(0, height/2, width, height/2);
  stroke(80, 255, 255);
  line(width/2, 0, width/2, height);
  translate(width/2, height/2);
  rotateY(PI/2);
  stroke(0, 255, 255);
  line(0, 0, -60, 0);
  popMatrix();
}
