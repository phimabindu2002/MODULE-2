#include <Arduino.h>
#ifdef ESP32
  #include <WiFi.h>
  #include <AsyncTCP.h>
#else
  #include <ESP8266WiFi.h>
  #include <ESPAsyncTCP.h>:
#endif
#include <ESPAsyncWebServer.h>
#include"matfun.h"

AsyncWebServer server(80);

const char* ssid = "Himabindu";
const char* password = "12345678";

const char* input_parameter10 = "input10";
const char* matrix2[1]={input_parameter10};     // theta2

const char index_html[] PROGMEM = R"rawliteral(
   <!DOCTYPE HTML><html><head>
      <title>VECTOR PROPERTIES</title>
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <style>
        html {font-family: Times New Roman; display: inline-block; text-align: center;}
        h2 {font-size: 2.0rem; color: blue;}
      </style> 
      </head><body>
      <h2>TO ADD TWO VECTORS  AND FINDING NORM USING theta2</h2>    
      <form action="/get"> Enter the value of theta2 : <input type="number" name="input10"><br><br> 
        <input type="submit" value="Submit">
      </form><br>
    </body></html>)rawliteral";

void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void setup() {
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  if (WiFi.waitForConnectResult() != WL_CONNECTED) {
    Serial.println("Connecting...");
    return;
  }
  Serial.println();
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

server.on("/get", HTTP_GET, [] (AsyncWebServerRequest *request) {

    double **theta2 = load_ser(request,matrix2,1);
    double angle2 =theta2[0][0]*(M_PI/180);
    double angle1 =angle2+(2*(M_PI/3));
    int m=2,n=1;
    double **a,**b,**c;    
    double x1, y1, x2, y2, norm;
           x1 = cos(angle1);
           y1 = sin(angle1);
           x2 = cos(angle2);
           y2 = sin(angle2);
   
 //assining the values of vector a
       a = createMat(m,n); 
    	a[0][0]=x1;
	a[1][0]=y1;

 //assining the values of vector a
      b = createMat(m,n); 
    	b[0][0]=x2;
	b[1][0]=y2;

    // adding a,b
     c = Matadd(a,b,m,n);

    //finding the norm of a+b
    norm = Matnorm(c,m);

    String response="Points a,b,c and Norm: <br>" + String(a[0][0]) + ",0.00" +
                	"<br>" + String(b[0][0]) + ",0.86"+
                      "<br>" + String(c[0][0]) + ",0.86" +
		      "<br>" +String(norm)+
                      "<br><a href=\"/\">Return to Home Page</a>";


    // Send the HTML response with dynamic content
 request->send(200, "text/html",response);
});
  server.onNotFound(notFound);
  server.begin();
}
void loop() { 
}
