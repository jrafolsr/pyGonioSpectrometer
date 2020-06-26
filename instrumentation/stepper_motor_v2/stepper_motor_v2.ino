/******************************************************************************
Modified version of the following demo:
Modified by Joan Ràfols, 11/06/2020.

Original code:
SparkFun Big Easy Driver Basic Demo
Toni Klopfenstein @ SparkFun Electronics
February 2015
https://github.com/sparkfun/Big_Easy_Driver

Simple demo sketch to demonstrate how 5 digital pins can drive a bipolar stepper motor,
using the Big Easy Driver (https://www.sparkfun.com/products/12859). Also shows the ability to change
microstep size, and direction of motor movement.

Development environment specifics:
Written in Arduino 1.6.0

This code is beerware; if you see me (or any other SparkFun employee) at the local, and you've found our code helpful, please buy us a round!
Distributed as-is; no warranty is given.

Example based off of demos by Brian Schmalz (designer of the Big Easy Driver).
http://www.schmalzhaus.com/EasyDriver/Examples/EasyDriverExamples.html
******************************************************************************/
//Declare pin functions on Arduino
#define gstep 2
#define gdir 3
#define MS1 4
#define MS2 5
#define MS3 6
#define EN  7
#define sdir 8 // Direction pin for the shutter motor
#define sstep 9 // Step pin for the shutter motor
#define sEN 10 // Enable the stepper of the shutter

#define DELAY 1

//Declare variables for functions
int rotation;
int x;
int y;
int state;
int steps;
int resolution;
int motor; // Which motor to move
int sopen =  1; // Variable that controls if the shutter is open or not



void setup() {
  pinMode(EN, OUTPUT);
  digitalWrite(EN, LOW);

  pinMode(sEN, OUTPUT);
  digitalWrite(sEN, HIGH);
  
  pinMode(gstep, OUTPUT);
  pinMode(gdir, OUTPUT);
  pinMode(MS1, OUTPUT);
  pinMode(MS2, OUTPUT);
  pinMode(MS3, OUTPUT);
  
  pinMode(sstep, OUTPUT);
  pinMode(sdir, OUTPUT);


  
  resetBEDPins(); //Set step, direction, microstep and enable pins to default states
  Serial.begin(9600); //Open Serial connection for debugging
  Serial.println("Motor ready");
  
  
}

//Main loop
void loop() {
  
  while(Serial.available() > 0){

      motor = Serial.parseInt(); //Read user input and trigger appropriate function
      
      if (motor == 0){
         steps = Serial.parseInt();
         
         if (Serial.read() == '\n'){
          MoveShutter(steps);
          Serial.println("Movement finished");
         }
      }
      else if (motor == 1){
        steps = Serial.parseInt();
        rotation = Serial.parseInt(); 
        resolution = Serial.parseInt();
      
/*      Serial.println(steps);
      Serial.println(rotation);
      Serial.println(resolution);*/
      
      
        if (Serial.read() == '\n') {
          if (rotation == 0)
          {
            resetBEDPins();// Reset the pins
            Serial.println("Movement finished");
          }
          else 
          {
          digitalWrite(EN, LOW); //Pull enable pin low to set FETs active and allow motor control
          SetStepResolution(resolution);
          delay(2); // It is said in the specs, to wait at least 1 ms before moving after the enable

            if (rotation == 1 || rotation == 2 )
            {
               MoveMotor(steps, rotation);
               Serial.println("Movement finished");
            }
            else
            {
              Serial.println("Invalid code string");
            }
          }
        }   
      }
    }
}

//Reset Big Easy Driver pins to default states
void resetBEDPins()
{
  digitalWrite(gstep, LOW);
  digitalWrite(gdir, LOW);
  digitalWrite(MS1, LOW);
  digitalWrite(MS2, LOW);
  digitalWrite(MS3, LOW);
  digitalWrite(EN, LOW);
  delay(1);
}
//Movement of the goniometer
void MoveMotor(int steps, int rotation){
  if (rotation == 1){
    digitalWrite(gdir, LOW); // Pull direction pin low to move "forward"
  }
  else if (rotation == 2){
    digitalWrite(gdir, HIGH); //Pull direction pin low to move in "reverse"
  }
  
  for(x= 0; x < steps; x++)  //Loop the stepping enough times for motion to be visible
  {
    digitalWrite(gstep,HIGH); //Trigger one step
    delay(DELAY);
    digitalWrite(gstep,LOW); //Pull step pin low so it can be triggered again
    delay(DELAY);

  }
}

// 1/16th microstep foward mode function
void SetStepResolution(int resolution)
{
  //Serial.println("Stepping at 1/16th microstep mode.");
  if (resolution == 1)
  {
    digitalWrite(MS1, LOW); //Pull MS1,MS2, and MS3  to set logic to 1/16th microstep resolution
    digitalWrite(MS2, LOW);
    digitalWrite(MS3, LOW);
  }
  else if (resolution == 2)
  {
    digitalWrite(MS1, HIGH); //Pull MS1,MS2, and MS3  to set logic to 1/2th microstep resolution
    digitalWrite(MS2, LOW);
    digitalWrite(MS3, LOW);
  }
  else if (resolution == 4)
  {
    digitalWrite(MS1, LOW); //Pull MS1,MS2, and MS3  to set logic to 1/4th microstep resolution
    digitalWrite(MS2, HIGH);
    digitalWrite(MS3, LOW);
  }
    else if (resolution == 8)
  {
    digitalWrite(MS1, HIGH); //Pull MS1,MS2, and MS3  to set logic to 1/8th microstep resolution
    digitalWrite(MS2, HIGH);
    digitalWrite(MS3, LOW);
  }  
  else if (resolution == 16)
  {
    digitalWrite(MS1, HIGH); //Pull MS1,MS2, and MS3 high to set logic to 1/16th microstep resolution
    digitalWrite(MS2, HIGH);
    digitalWrite(MS3, HIGH);
  }
  else
  {
    digitalWrite(MS1, LOW); //Pull MS1,MS2, and MS3  to set logic to 1/16th microstep resolution
    digitalWrite(MS2, LOW);
    digitalWrite(MS3, LOW);
  }
}

void MoveShutter(int steps){
  // Enabling the shutter motor
  digitalWrite(sEN, LOW);
  delay(2);
  
//Serial.println("Moving the shutter according to the global "opens" ");
  if (sopen == 0){
    digitalWrite(sdir, LOW); //Pull direction pin low to move "forward"
    sopen = 1; // Change the variable so next time in does the reverse
  }
  else if (sopen == 1){
    digitalWrite(sdir, HIGH); //Pull direction pin low to move in "reverse"
    sopen = 0;
  }
  
  for(x= 0; x < steps; x++)  // Move 90 degrees in the define direction, assuming a 1.8deg step
  {
    digitalWrite(sstep,HIGH); //Trigger one step forward
    delay(DELAY);
    digitalWrite(sstep,LOW); //Pull step pin low so it can be triggered again
    delay(DELAY);
    }
  delay(2000);
  // Disabling the shutter motor
  digitalWrite(sEN, HIGH);
  delay(2);
  
}
