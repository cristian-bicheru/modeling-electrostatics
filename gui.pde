/* =========================================================
 * ====                   WARNING                        ===
 * =========================================================
 * The code in this tab has been generated from the GUI form
 * designer and care should be taken when editing this file.
 * Only add/edit code inside the event handlers i.e. only
 * use lines between the matching comment tags. e.g.

 void myBtnEvents(GButton button) { //_CODE_:button1:12356:
     // It is safe to enter your event code here  
 } //_CODE_:button1:12356:
 
 * Do not rename this tab!
 * =========================================================
 */

public void button2_click1(GButton source, GEvent event) { //_CODE_:button2:309345:
  println("button2 - GButton >> GEvent." + event + " @ " + millis());
  ocx = ox;
  ocy = oy;
  ocz = oz;
  dex = 0;
  dey = 0;
  dez = 0;
  slider1.setValue(0.5);
  slider2.setValue(0.5);
  slider3.setValue(0.5);
} //_CODE_:button2:309345:

public void checkbox3_clicked1(GCheckbox source, GEvent event) { //_CODE_:checkbox3:661901:
  println("checkbox3 - GCheckbox >> GEvent." + event + " @ " + millis());
  showField = (showField + 1) % 2;
  showMagneticField = (showMagneticField + 1) % 2;
} //_CODE_:checkbox3:661901:

public void slider1_change1(GSlider source, GEvent event) { //_CODE_:slider1:735542:
  println("slider1 - GSlider >> GEvent." + event + " @ " + millis());
  dex = (source.getValueF()-0.5)*720;
} //_CODE_:slider1:735542:

public void slider2_change1(GSlider source, GEvent event) { //_CODE_:slider2:395567:
  println("slider2 - GSlider >> GEvent." + event + " @ " + millis());
  dey = (source.getValueF()-0.5)*720;
} //_CODE_:slider2:395567:

public void slider3_change1(GSlider source, GEvent event) { //_CODE_:slider3:966413:
  println("slider3 - GSlider >> GEvent." + event + " @ " + millis());
  dez = (source.getValueF()-0.5)*720;
} //_CODE_:slider3:966413:



// Create all the GUI controls. 
// autogenerated do not edit
public void createGUI(){
  G4P.messagesEnabled(false);
  G4P.setGlobalColorScheme(GCScheme.BLUE_SCHEME);
  G4P.setMouseOverEnabled(false);
  surface.setTitle("Sketch Window");
  label1 = new GLabel(this, 70, 40, 120, 20);
  label1.setTextAlign(GAlign.CENTER, GAlign.MIDDLE);
  label1.setText("Simulation Controls");
  label1.setLocalColorScheme(GCScheme.SCHEME_15);
  label1.setOpaque(false);
  button2 = new GButton(this, 90, 80, 80, 20);
  button2.setText("Reset Graph");
  button2.setLocalColorScheme(GCScheme.SCHEME_15);
  button2.addEventHandler(this, "button2_click1");
  checkbox3 = new GCheckbox(this, 60, 110, 140, 20);
  checkbox3.setIconAlign(GAlign.LEFT, GAlign.MIDDLE);
  checkbox3.setText("Show Inputted Field");
  checkbox3.setLocalColorScheme(GCScheme.SCHEME_15);
  checkbox3.setOpaque(false);
  checkbox3.addEventHandler(this, "checkbox3_clicked1");
  slider1 = new GSlider(this, 80, 130, 100, 40, 10.0);
  slider1.setLimits(0.5, 0.0, 1.0);
  slider1.setNumberFormat(G4P.DECIMAL, 2);
  slider1.setOpaque(false);
  slider1.addEventHandler(this, "slider1_change1");
  slider2 = new GSlider(this, 80, 160, 100, 40, 10.0);
  slider2.setLimits(0.5, 0.0, 1.0);
  slider2.setNumberFormat(G4P.DECIMAL, 2);
  slider2.setOpaque(false);
  slider2.addEventHandler(this, "slider2_change1");
  slider3 = new GSlider(this, 80, 190, 100, 40, 10.0);
  slider3.setLimits(0.5, 0.0, 1.0);
  slider3.setNumberFormat(G4P.DECIMAL, 2);
  slider3.setOpaque(false);
  slider3.addEventHandler(this, "slider3_change1");
}

// Variable declarations 
// autogenerated do not edit
GLabel label1; 
GButton button2; 
GCheckbox checkbox3; 
GSlider slider1; 
GSlider slider2; 
GSlider slider3; 
