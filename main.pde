import g4p_controls.*;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

int width = 1600;
int height = 1000;

int planeMode = 0;
int colorMode = 0;
float scaleFactor = 0.01; // m/pixels


//plane variables
float px, py, pz, pt, pr = 0;

//camera stuff
float dez = 0;
float dex = -40;
float dey = -50;
float ox = width/2+150;
float oy = height/2;
float oz = 0;
float ocx = ox;
float ocy = oy;
float ocz = oz;

//vector functions
String fi, fj, fk;

//array of all wires
wire[] wires = {};

//vector field
float[][] vectorField = {};
int vectorDetail = 10;
int showField = 0;
float vectorScale = 1;

//net magnetic field
float[][] magneticField = {};
int showMagneticField = 1;
float magneticVectorScale = 1;
double dt = 0.1;
double tension = 1;
double[][] splineMatrix = {{0, 1, 0, 0}, {-tension, 0, tension, 0}, {2*tension, tension-3, 3-2*tension, -tension}, {-tension, 2-tension, tension-2, tension}};

//field lines
float[][] startingPoints = {};
float[][] lines = {};

//create eval engine
ScriptEngineManager mgr = new ScriptEngineManager();
ScriptEngine engine = mgr.getEngineByName("javascript");

void setupEvalEngine() {
  try {
    engine.eval("abs = Math.abs; acos = Math.acos; acosh = Math.acosh; asin = Math.asin; asinh = Math.asinh; atan = Math.atan; atan2 = Math.atan2; atanh = Math.atanh; cbrt = Math.cbrt; ceil = Math.ceil; clz32 = Math.clz32; cos = Math.cos; cosh = Math.cosh; exp = Math.exp; expm1 = Math.expm1; floor = Math.floor; fround = Math.fround; hypot = Math.hypot; imul = Math.imul; log = Math.log; log1p = Math.log1p; log2 = Math.log2; log10 = Math.log10; maax = Math.maax; min = Math.min; pow = Math.pow; random = Math.random; round = Math.round; sign = Math.sign; sin = Math.sin; sinh = Math.sinh; sqrt = Math.sqrt; tan = Math.tan; tanh = Math.tanh; trunc = Math.trunc;");
  } catch (ScriptException e) {
     e.printStackTrace();
  }
}

//vector function evaluators
float i(float x, float y, float z) {
  String evalString;
  try {
    evalString = fi.replace("x", str(x)).replace("y", str(y)).replace("z", str(z));
    engine.eval("fi = "+evalString+';');
    return (float)(double)engine.get("fi");
  } catch (Exception e) {
    try {
      return (float)(int)engine.get("fi");
    } catch (Exception f) {
    }
    return 0;
  }
}

float j(float x, float y, float z) {
  String evalString;
  try {
    evalString = fj.replace("x", str(x)).replace("y", str(y)).replace("z", str(z));
    engine.eval("fj = "+evalString+';');
    return (float)(double)engine.get("fj");
  } catch (Exception e) {
    try {
      return (float)(int)engine.get("fj");
    } catch (Exception f) {
    }
    return 0;
  }
}

float k(float x, float y, float z) {
  String evalString;
  try {
    evalString = fk.replace("x", str(x)).replace("y", str(y)).replace("z", str(z));
    engine.eval("fk = "+evalString+';');
    return (float)(double)engine.get("fk");
  } catch (Exception e) {
    try {
      return (float)(int)engine.get("fk");
    } catch (Exception f) {
    }
    return 0;
  }
}
//

int test = 0;

void setup() {
  size(1600, 1000, P3D);
  createGUI();
  importScene();
  setupEvalEngine();
  calculateVectorField(-200, 200, -200, 200, -200, 200);
  //calculateNetMagneticField(-150, 150, -150, 150, -100, 100);
  calculateLines();
  //dump();
}

void dump() {
  String[] strcoords = {};
    for (int itr = 0; itr < lines.length; itr++) {
      String[] vctrstr = {Float.toString(lines[itr][0]), Float.toString(lines[itr][1]), Float.toString(lines[itr][2]), Float.toString(lines[itr][3]), Float.toString(lines[itr][4]), Float.toString(lines[itr][5])};
      strcoords = concat(strcoords, vctrstr);
    }
    saveStrings("dump.txt", strcoords);
}

void draw() {
  background(102);
  pushMatrix();
  translate(ocx, ocy, ocz);
  rotateX(dey/100);
  rotateY(dex/100);
  rotateZ(dez/100);
  drawGraph(-240, 240, -240, 240, -240, 240);
  if (showField == 1) {
    drawVectorField();
  }
  if (showMagneticField == 2) {
    drawMagneticField();
  }
  drawLines();
  for (int i = 0; i < wires.length; i++) {
    wires[i].draw();
  }
  popMatrix();
}

double[] crsdot(double t) {
  double [] result = {};
  double [] a = {1, t, t*t, t*t*t};
  double temp = 0;
  for (int c = 0; c < 4; c++) {
    temp = 0;
    for (int i = 0; i < 4; i++) {
      temp += a[i]*splineMatrix[i][c];
    }
    result = (double[])append(result, temp);
  }
  return result;
}

PVector integrateCatmullRomSpline(float cx, float cy, float cz, float x1, float y1, float z1, float x2, float y2, float z2, float cx2, float cy2, float cz2, float ppx, float ppy, float ppz) {
  PVector result = new PVector(0, 0, 0);
  double[] dotResult = {};
  double xi, yi, zi, xf, yf, zf;
  float mag;
  PVector ds, r, tVector;
  for (double t = 0; t < 1-dt; t += dt) {/**
    dotResult = crsdot(t);
    xi = cx*dotResult[0]+x1*dotResult[1]+x2*dotResult[2]+cx2*dotResult[3];
    yi = cy*dotResult[0]+y1*dotResult[1]+y2*dotResult[2]+cy2*dotResult[3];
    zi = cz*dotResult[0]+z1*dotResult[1]+z2*dotResult[2]+cz2*dotResult[3];
    dotResult = crsdot(t+dt);
    xf = cx*dotResult[0]+x1*dotResult[1]+x2*dotResult[2]+cx2*dotResult[3];
    yf = cy*dotResult[0]+y1*dotResult[1]+y2*dotResult[2]+cy2*dotResult[3];
    zf = cz*dotResult[0]+z1*dotResult[1]+z2*dotResult[2]+cz2*dotResult[3]; **/
    xi = curvePoint(x1, cx, cx2, x2, (float) t);
    yi = curvePoint(y1, cy, cy2, y2, (float) t);
    zi = curvePoint(z1, cz, cz2, z2, (float) t);
    xf = curvePoint(x1, cx, cx2, x2, (float) (t+dt));
    yf = curvePoint(y1, cy, cy2, y2, (float) (t+dt));
    zf = curvePoint(z1, cz, cz2, z2, (float) (t+dt));
    
    ds = new PVector((float) (xf-xi), (float) (yf-yi), (float) (zf-zi));
    r = new PVector((float) (ppx-xi), (float) (ppy-yi), (float) (ppz-zi));
    mag = r.magSq();
    r.normalize();
    tVector = ds.cross(r);
    tVector.div(mag);
    result = result.add(tVector); 
  }
  
  return result;
}

void calculateNetMagneticField(float x1, float x2, float y1, float y2, float z1, float z2) {
  //only for large scale, needs method for when min < 16
  int delta, len, idx = 0;
  float sx, sy, sz, ti, tj, tk, max = 0;
  PVector fieldVector;
  PVector tVector;
  float[][] wireCoords;
  delta = 40;
  for (float x = ceil(x1); x <= ceil(x2); x+= delta) {
    for (float y = ceil(y1); y <= ceil(y2); y+= delta) {
      for (float z = ceil(z1); z <= ceil(z2); z+= delta) {
        sx = x*scaleFactor;
        sy = y*scaleFactor;
        sz = z*scaleFactor;
        ti = i(sx, sy, sz);
        tj = j(sx, sy, sz);
        tk = k(sx, sy, sz);
        fieldVector = new PVector(0,0,0);
        for (int i = 0; i < wires.length; i++) {
          wireCoords = wires[i].getCoords();
          len = wireCoords.length;
          if (len > 4) { // broken, calculateLines has the correct code for this
            for (int j = 0; j < (len-(len % 4)); j += 4) {
              tVector = integrateCatmullRomSpline(wireCoords[j][0]*scaleFactor, wireCoords[j][1]*scaleFactor, wireCoords[j][2]*scaleFactor, wireCoords[j+1][0]*scaleFactor, wireCoords[j+1][1]*scaleFactor, wireCoords[j+1][2]*scaleFactor, wireCoords[j+2][0]*scaleFactor, wireCoords[j+2][1]*scaleFactor, wireCoords[j+2][2]*scaleFactor, wireCoords[j+3][0]*scaleFactor, wireCoords[j+3][1]*scaleFactor, wireCoords[j+3][2]*scaleFactor, sx, sy, sz);
              tVector.mult(wires[i].getCurrent()/(4*PI));
              fieldVector = fieldVector.add(tVector);
            }
          }
          tVector = integrateCatmullRomSpline(wireCoords[len-4][0]*scaleFactor, wireCoords[len-4][1]*scaleFactor, wireCoords[len-4][2]*scaleFactor, wireCoords[len-3][0]*scaleFactor, wireCoords[len-3][1]*scaleFactor, wireCoords[len-3][2]*scaleFactor, wireCoords[len-2][0]*scaleFactor, wireCoords[len-2][1]*scaleFactor, wireCoords[len-2][2]*scaleFactor, wireCoords[len-1][0]*scaleFactor, wireCoords[len-1][1]*scaleFactor, wireCoords[len-1][2]*scaleFactor, sx, sy, sz);
          tVector.mult(wires[i].getCurrent()/(4*PI));
          fieldVector = fieldVector.add(tVector);
        }
        float[] vals = {max, abs(fieldVector.x+ti), abs(fieldVector.y+tj), abs(fieldVector.z+tk)};
        max = max(vals);
        float[] vector = {x, y, z, fieldVector.x+ti, fieldVector.y+tj, fieldVector.z+tk};
        magneticField = (float[][])append(magneticField, vector);
      }
    }
  }
  magneticVectorScale = 100/max;
  println("scale");
  println(magneticVectorScale);
}

void calculateLines() {
  float sx, sy, sz, ti, tj, tk, x, y, z, initx, inity, initz, sdist = 21;
  float dist = 0;
  PVector fieldVector, tVector;
  int len;
  float[][] wireCoords;
  for (int l = 0; l<startingPoints.length; l++) {
    float[][] line = {};
    x = startingPoints[l][0];
    y = startingPoints[l][1];
    z = startingPoints[l][2];
    sx = x*scaleFactor;
    sy = y*scaleFactor;
    sz = z*scaleFactor;
    initx = x;
    inity = y;
    initz = z;
    dist = 0;
    sdist = 0;
    while (dist<4000 && (sdist>200 || dist<100)) {
      ti = i(sx, sy, sz);
      tj = j(sx, sy, sz);
      tk = k(sx, sy, sz);
      fieldVector = new PVector(0,0,0);
      for (int i = 0; i < wires.length; i++) {
        wireCoords = wires[i].getCoords();
        len = wireCoords.length;
        for (int j = 0; j < len-3; j++) {
          tVector = integrateCatmullRomSpline(wireCoords[j][0]*scaleFactor, wireCoords[j][1]*scaleFactor, wireCoords[j][2]*scaleFactor, wireCoords[j+1][0]*scaleFactor, wireCoords[j+1][1]*scaleFactor, wireCoords[j+1][2]*scaleFactor, wireCoords[j+2][0]*scaleFactor, wireCoords[j+2][1]*scaleFactor, wireCoords[j+2][2]*scaleFactor, wireCoords[j+3][0]*scaleFactor, wireCoords[j+3][1]*scaleFactor, wireCoords[j+3][2]*scaleFactor, sx, sy, sz);
          tVector.mult(wires[i].getCurrent()/(4*PI));
          fieldVector = fieldVector.add(tVector);
        }
      }
      fieldVector.normalize();
      dist += 1;
      float[] vect = {x, y, z, fieldVector.x+ti, fieldVector.y+tj, fieldVector.z+tk};
      lines = (float[][])append(lines, vect);
      x += (fieldVector.x+ti);
      y += (fieldVector.y+tj);
      z += (fieldVector.z+tk);
      sx = x*scaleFactor;
      sy = y*scaleFactor;
      sz = z*scaleFactor;
      sdist = pow((x-initx), 2)+pow((y-inity), 2)+pow((z-initz), 2);
    }
  }
}

void drawLines() {
  int mod = 0;
  fill(#336EFF);
  stroke(#336EFF); 
  for (int i = 0; i < lines.length; i++) {
    if (mod%100 == 0) {
       drawVector(lines[i], 1);
       mod = 1;
    } else {
       line(-lines[i][0], -lines[i][1], -lines[i][2], -lines[i][0]-lines[i][3], -lines[i][1]-lines[i][4], -lines[i][2]-lines[i][5]);
       mod++;
    }
  }
}

void drawMagneticField() {
  fill(#f93e2c);
  stroke(#f93e2c);
  strokeWeight(3);
  for (int i = 0; i < magneticField.length; i++) {
    drawVector(magneticField[i], magneticVectorScale);
  }
}

void mouseDragged() {
  if (mouseButton == CENTER) {
    ocx += mouseX-pmouseX;
    ocy += mouseY-pmouseY;
  }
}

void mouseWheel(MouseEvent event) {
  float e = event.getCount();
  ocz += -e*20;
}

void calculateVectorField(float x1, float x2, float y1, float y2, float z1, float z2) {
  //only for large scale, needs method for when min < 16
  int delta, idx = 0;
  float sx, sy, sz, ti, tj, tk, max = 0;
  delta = floor(min(z2-z1, x2-x1)/4);
  for (float x = ceil(x1); x <= ceil(x2); x+= delta) {
    for (float y = ceil(y1); y <= ceil(y2); y+= delta) {
      for (float z = ceil(z1); z <= ceil(z2); z+= delta) {
        sx = x*scaleFactor;
        sy = y*scaleFactor;
        sz = z*scaleFactor;
        ti = i(sx, sy, sz);
        tj = j(sx, sy, sz);
        tk = k(sx, sy, sz);
        float[] vals = {max, ti, tj, tk};
        max = max(vals);
        float[] vector = {x, y, z, ti, tj, tk};
        vectorField = (float[][])append(vectorField, vector);
      }
    }
  }
  if (max != 0) {
    vectorScale = 20/max;
  } else {
    vectorScale = 0;
  }
}

void drawCone(float r, float h) {
  circle(0, 0, 2*r);
  float theta = 0;
  float dtheta = 2*PI/vectorDetail;
  beginShape(TRIANGLES);
  for (int idx = 0; idx < vectorDetail; idx++) {
    vertex(0, 0, h);
    vertex(r*sin(theta), r*cos(theta), 0);
    vertex(r*sin(theta+dtheta), r*cos(theta+dtheta), 0);
    theta += dtheta;
  }
  endShape();
}

float closestTo0(float a, float b) {
  if (abs(a) < abs(b)) {
    return a;
  }
  return b;
}

void drawVector(float[] vector, float scale) {
  float x = vector[0];
  float y = vector[1];
  float z = vector[2];
  float i = vector[3];
  float j = vector[4];
  float k = vector[5];
  float lx = i*scale;
  float ly = j*scale;
  float lz = k*scale;
  if (pow(lx, 2)+ pow(ly, 2)+ pow(lz, 2) > 400) {
    PVector tVector;
    tVector = new PVector(lx, ly, lz);
    tVector.normalize();
    tVector.mult(20);
    lx = tVector.x;
    ly = tVector.y;
    lz = tVector.z;
  }
  float dx = sqrt(pow(i, 2)+pow(k, 2));
  pushMatrix();
  translate(-x, -y, -z);
  line(0, 0, 0, -lx, -ly, -lz);
  translate(-lx, -ly, -lz);
  rotateY(atan2(-i, -k));
  rotateX(atan2(j, dx));
  drawCone(3, 5);
  popMatrix();
}

void drawVectorField() {
  fill(#f93e2c);
  stroke(#f93e2c);
  strokeWeight(3);
  for (int i = 0; i < vectorField.length; i++) {
    drawVector(vectorField[i], vectorScale);
  }
}

void drawGraph(float x1, float x2, float y1, float y2, float z1, float z2) {
  //only for large scale, needs method for when min < 16
  int delta, x, z;
  fill(255);
  stroke(#B52B2B);
  strokeWeight(1);
  line(x1, 0, 0, x2, 0, 0);
  line(0, y1, 0, 0, y2, 0);
  line(0, 0, z1, 0, 0, z2);
  delta = floor(min(z2-z1, x2-x1)/16);
  if (delta == 0) {
    println("no");
  }
  textSize(delta/2);
  textAlign(CENTER);
  x = ceil(x1);
  while (true) {
    line(x, 0, z1, x, 0, z2);
    x += delta;
    if (x > x2) {
      break;
    }
  }
  pushMatrix();
  rotateY(PI/2);
  x = ceil(x1);
  while (true) {
    text(str(x), 0, 0, x);
    x += delta;
    if (x > x2) {
      break;
    }
  }
  popMatrix();
  z = ceil(z1);
  while (true) {
    line(x1, 0, z, x2, 0, z);
    text(str(z), 0, 0, z);
    z += delta;
    if (z > z2) {
      break;
    }
  }
}

//kind of uneeded, the -y is though
float[] transformCoords(float[] point) {
  float[] ret = {point[0], -point[1], point[2]};
  return ret;
}

class wire {
  float[][] coords = {};
  float current;
  
  wire(float[][] tcoords, float tcurrent) {
    coords = (float[][])concat(coords, tcoords);
    current = tcurrent;
  }
  
  void draw() {
    float[] p1, p2, p3, p4;
    stroke(#FFD700);
    strokeWeight(4);
    strokeCap(PROJECT);
    noFill();
    for (int i = 0; i < (coords.length-3); i++) {
      p1 = transformCoords(coords[i]);
      p2 = transformCoords(coords[i+1]);
      p3 = transformCoords(coords[i+2]);
      p4 = transformCoords(coords[i+3]);
      curve(p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2], p4[0], p4[1], p4[2]);
    }
  }
  
  float[][] getCoords() {
    return coords;
  }
  
  float getCurrent() {
    return current;
  }
}

void importScene() {
  String cLine;
  String[] sceneDataRaw = loadStrings("scenetransformer.mgs");
  for (int i = 0; i < sceneDataRaw.length; i++) {
    cLine = sceneDataRaw[i];
    if (cLine.equals("field {")) {
      fi = sceneDataRaw[i+1];
      fj = sceneDataRaw[i+2];
      fk = sceneDataRaw[i+3];
      i += 3;
    } else if (cLine.equals("wire {")) {
      float[][] coords = {};
      float current;
      i++;
      for (int j = i; j < sceneDataRaw.length; j++) {
        cLine = sceneDataRaw[i];
        if (cLine.indexOf(",") != -1) {
          coords = (float[][])append(coords, float(split(cLine, ", ")));
        } else {
          break;
        }
        i++;
      }
      current = float(sceneDataRaw[i]);
      wires = (wire[])append(wires, new wire(coords, current));
      i++;
    } else if (cLine.equals("starts {")) {
      i++;
      for (int j = i; j < sceneDataRaw.length; j++) {
        cLine = sceneDataRaw[i];
        if (cLine.indexOf(",") != -1) {
          startingPoints = (float[][])append(startingPoints, float(split(cLine, ", ")));
        } else {
          break;
        }
        i++;
      }   
    }
  }
}
