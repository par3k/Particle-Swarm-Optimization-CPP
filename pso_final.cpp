/*
    ########    #######   ########
    ##     ## ##      ## ##      ##
    ##     ## ##         ##      ##
    ########   #######   ##      ##
    ##                ## ##      ##
    ##         ##     ## ##      ##
    ##          #######   ########  

Final project of the C++ course
18309002 - PARK HOIJAI

Implement PSO with C++
    Higher scores if you implement advanced PSO
Solve several optimization problems
    At most 15 problems
    Solve more pronlems, get higher scores
Required materials
    Document - analysis of your codes, solved problems, Figures of final results, Comments about of PSO
    PPT files
    Source code
Deadline
    18:00 on 11/6
*/

#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#define rand_01 ((float)rand() / (float)RAND_MAX)

// set dimension and number of particles
const int numofdims = 10;
const int numofparticles = 20;

using namespace std;

// Sphere function Equation
void F1(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
    float t1;
	memset(fitnesses, 0, sizeof (float) * numofparticles); // create space fitnesses[i] array
    for(int i = 0; i < numofparticles; i++){
        for(int j = 0; j < numofdims; j++){
        	t1 = X[i][j];
        	t1 *= t1;
            fitnesses[i] += t1; 
        }
    }
    p= fitnesses[numofparticles-1];
}

// Schwefel's 2.22 function Equation
void F2(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float t1,t2;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		t1 = abs(X[i][j]);
        		t2 *= abs(X[i][j]);
        		fitnesses[i] += t1 ;//+ t2;
		}
	}
	    p= fitnesses[numofparticles-1];
}

// Schwefel's 1.20 function Equation
void F3(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float t1;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		t1 = X[i][j-1];
        		t1 += t1;
        		t1 *= t1;
        		fitnesses[i] += t1;
		}
	}
	p= fitnesses[numofparticles-1];
}

// Schwefel's 2.21 function Equation
void F4(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float t1;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		t1 = abs(X[i][j]);
        		t1 *= t1;
        		fitnesses[i] = t1;
		}
	}
	p= fitnesses[numofparticles-1];
}

// Rosenbrock function Equation
void F5(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float x1, x2, t1, t2;
	memset(fitnesses, 0, sizeof (float) * numofparticles);    
    for(int i = 0; i < numofparticles; i++)
        for(int j = 0; j < numofdims - 1; j++){
            x1 = X[i][j];
            x2 = X[i][j+1];
            t1 = (x2 - x1 * x1);
            t1 *= t1;
            t1 *= 100;
            t2 = x1 - 1;
            t2 *= t2;
            fitnesses[i] = t1 + t2;
        }
        p= fitnesses[numofparticles-1];
}

// Step function Equation
void F6(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float t1;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		t1 = X[i][j] + 0.5;
        		t1 *= t1;
        		fitnesses[i] += t1;
		}
	}
	p= fitnesses[numofparticles-1];
}

// Quartic noise function Equation
void F7(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float x1,t1,t2;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		x1 = X[i][j];
        		x1 *= x1;
        		x1 *= x1;
        		x1 *= j;
        		t1 += x1;
        		t2 = rand() % 2;
        		fitnesses[i] = t1 + t2; 
        		
		}
	}
}

// Schwefel function Equation
void F8(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		fitnesses[i] = ((- X[i][j])* sin (sqrt(abs( X[i][j] ))));
		}
	}
}

// Rastrigin function Equation
void F9(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	double pi = acos(-1);
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		fitnesses[i] += (pow(X[i][j],2) - 10 * cos (2 * pi * X[i][j]) + 10);
		}
	}
}

// Ackley function Equation
void F10(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float x1,x2,t1,t2;
	double pi = acos(-1);
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		x1 = X[i][j];
        		x1 *= x1;
        		t1 = x1 / j;
        		t2 = cos(2 * pi * X[i][j]) / j;
        		fitnesses[i] = -20 * exp(-0.2 * sqrt(t1)) - exp(t2) + 20 + exp(1);
		}
	}
}

// Griewank function Equation
void F11(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float t1;
	double pi = acos(-1);
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims; j++){
        		t1 = 1;
        		fitnesses[i] = pow(X[i][j],2) / 4000 + 1;
        		for(int k = 0; k < numofdims; k++){
        			t1 = t1 * cos(X[i][k] / sqrt(k+1));
				}
				fitnesses[i] = fitnesses[i] - t1;
		}
	}
}

void F12(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float y1,y2,y3,y4,t1,t2,t3,t4;
	float x1,k1,k2;
	float z;
	double pi = acos(-1);
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < numofdims-1; j++){
        		y1 = X[i][1];
        		t1 = pow(sin(pi * y1),2);
        		t1 *= 10;
				y2 = X[i][j];
        		t2 = pow((y2 - 1),2);
				y3 = X[i][j+1];
				t3 = pow(sin ( pi * y3),2);
				t3 *= 10;
				y4 = X[i][numofdims];
				t4 = pow((y4 - 1),2);
				fitnesses[i] = t2 * ( 1 + t3 + t4);
		
				for (int h = 0; h < numofdims; h++){
					if(X[i][h] > 10){
						z = 100 * (X[i][h]-10) * (X[i][h]-10) * (X[i][h]-10) * (X[i][h]-10);
					} else if(-10 <= X[i][h] <= 10){
						z = 0;
					} else {
						z = 100 * (- X[i][h] - 10) * (- X[i][h] - 10) * (- X[i][h] - 10) * (- X[i][h] - 10);
					}
				}	
				fitnesses[i] = ((pi / numofdims)* (t1 * fitnesses[i]))+ z;			
			}	
				
	}
}

void F13(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float x1,x2,t1,t2,t3,t4,t5,t6;
	float z;
	double pi = acos(-1);
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        		x1 = X[i][1];
        		t1 = sin ( 3 * pi * x1);
        		t1 *= t1;	    	
        	for(int j = 0; j < numofdims; j++){
        		x2 = X[i][j];
        		t2 = x2 - 1;
        		t2 *= t2;
        		t3 = sin ( 3 * pi * x2 + 1);
        		t3 *= t3;
        		
        		fitnesses[i] += t2 * (1 + t3);
        		
        		t4 = x2 - 1;
        		t4 *= t4;
        		t5 = sin ( 2 * pi * x2);
        		t5 *= t5;
        		t6 = t4 * t5;
        		
				for (int h = 0; h < numofdims; h++){
					if(X[i][h] > 5){
						z = 100 * (X[i][h]-5) * (X[i][h]-5) * (X[i][h]-5) * (X[i][h]-5);
					} else if(-5 <= X[i][h] <= 5){
						z = 0;
					} else {
						z = 100 * (- X[i][h] - 5) * (- X[i][h] - 5) * (- X[i][h] - 5) * (- X[i][h] - 5);
					}
				}	        		
		fitnesses[i] = 0.1 * (t1 + fitnesses[i] + t6) + z;
		}

	}
}

void F14(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float x1,a1,t1;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < 25; j++){
        		for(int k = 0; k < 2; k++){
        			x1 = X[i][k];
        			a1 = X[k][j];
        			t1 = x1 - a1;
        			t1 = t1 * t1 * t1 * t1 * t1 * t1;
        			fitnesses[i] += t1;
        			fitnesses[i] = j + fitnesses[i];
        			fitnesses[i] = 1 / fitnesses[i];
        			}
        	fitnesses[i] += fitnesses[i];
		}
		fitnesses[i] = 1/500 + fitnesses[i];
		fitnesses[i] = 1 / fitnesses[i];
	}
}

void F15(float X[numofparticles][numofdims], float fitnesses[numofparticles], float &p){ 
	float a1,b1,t1,t2,t3,t4;
	memset(fitnesses, 0, sizeof (float) * numofparticles);
	    for(int i = 0; i < numofparticles; i++){
        	for(int j = 0; j < 11; j++){
        		a1 = X[1][j];
        		b1 = X[2][j];
        		t1 = X[i][j] * (b1 * b1 + b1 * X[i][2]);
        		t2 = (b1 * b1 + b1 * X[i][3] + X[i][4]);
        		t3 = t1 / t2;
        		t4 = a1 - t3;
        		t4 *= t4;
        		fitnesses[i] = t4;
		}
	}

}


//Particle Swarm Optimization **************
void PSO(int numofiterations, float c1, float c2, 
              float Xmin[numofdims], float Xmax[numofdims], float initialpop[numofparticles][numofdims],
                float bests[], float *gbestfit, float gbest[numofdims], float* p, int select){
    
    // STEP 1. value declaration
	float V[numofparticles][numofdims] = {0};
    float X[numofparticles][numofdims];
    float Vmax[numofdims];
    float Vmin[numofdims];
    float pbests[numofparticles][numofdims];
    float pbestfits[numofparticles];
    float fitnesses[numofparticles];
    float w;
    float minfit;
    int   minfitidx;
    
    //copies initial random values array to initial positions of the particles
    memcpy(X, initialpop, sizeof(float) * numofparticles * numofdims);
    
    	switch(select){
		case 1:
			F1(X, fitnesses, *p);
			break;
		case 2:
			for (int i = 0; i < numofparticles; i++)
				{
				Xmax[i] = 10;
				Xmin[i] = -10;
				}
			F2(X, fitnesses, *p);
			break;
		case 3:
			F3(X, fitnesses, *p);
			break;
		case 4:
			F4(X, fitnesses, *p);
			break;
		case 5:
			F5(X, fitnesses, *p);
			break;
		case 6:
			for (int i = 0; i < numofparticles; i++)
				{
				Xmax[i] = 1.28;
				Xmin[i] = -1.28;
				}
			F6(X, fitnesses, *p);
			break;
		case 7:
			F7(X, fitnesses, *p);
			break;
		case 8:
			for (int i = 0; i < numofparticles; i++)
				{
				Xmax[i] = 500;
				Xmin[i] = -500;
				}
			F8(X, fitnesses, *p);
			break;
		case 9:
			F9(X, fitnesses, *p);
			break;
		case 10:
			F10(X, fitnesses, *p);
			break;
		case 11:
			F11(X, fitnesses, *p);
			break;
		case 12:
			F12(X, fitnesses, *p);
			break;
		case 13:
			F13(X, fitnesses, *p);
			break;
		case 14:
			F14(X, fitnesses, *p);
			break;
		case 15:
			F15(X, fitnesses, *p);
			break;
		default:
			;
		}
    
    
    //F1(X, fitnesses); // Just change this F*(X, fitnesses) part, if you want to use other questions.
    minfit = *min_element(fitnesses, fitnesses + numofparticles);
    minfitidx = min_element(fitnesses, fitnesses + numofparticles) - fitnesses;
    *gbestfit = minfit; // gbestfit initialize
    memcpy(gbest, X[minfitidx], sizeof(float) * numofdims); // end copy
    
	// STEP 2. Setting speed limit
    for(int i = 0; i < numofdims; i++){ 
        Vmax[i] = 0.2 * (Xmax[i] - Xmin[i]);
        Vmin[i] = -Vmax[i];
    }

	// STEP 3. Calculate the minimum value of individual history
	// loop iterates 1000 times through the pso algorithm
    for(int t = 0; t < 1000; t++){ 
        w = 0.9 - 0.7 * t / numofiterations;

		// STEP 4. set pbest - The best position of the particle
        for(int i = 0; i < numofparticles; i++){
            if(fitnesses[i] < pbestfits[i]){
                pbestfits[i] = fitnesses[i];  
                memcpy(pbests[i], X[i], sizeof(float) * numofdims);
            }
        }
       
	    // STEP 5. Particle update rule
        for(int i = 0; i < numofparticles; i++){
            for(int j = 0; j < numofdims; j++){
            	// velocity
                V[i][j] = min(max((w * V[i][j] + rand_01 * c1 * (pbests[i][j] - X[i][j]) + rand_01 * c2 * (gbest[j] - X[i][j])), Vmin[j]), Vmax[j]);
                // position : p = p + v
                X[i][j] = min(max((X[i][j] + V[i][j]), Xmin[j]), Xmax[j]);
            }
        }

    	switch(select){
		case 1:
			F1(X, fitnesses, *p);
			break;
		case 2:
			for (int i = 0; i < numofparticles; i++)
				{
				Xmax[i] = 10;
				Xmin[i] = -10;
				}
			F2(X, fitnesses, *p);
			break;
		case 3:
			F3(X, fitnesses, *p);
			break;
		case 4:
			F4(X, fitnesses, *p);
			break;
		case 5:
			F5(X, fitnesses, *p);
			break;
		case 6:
			for (int i = 0; i < numofparticles; i++)
				{
				Xmax[i] = 1.28;
				Xmin[i] = -1.28;
				}
			F6(X, fitnesses, *p);
			break;
		case 7:
			F7(X, fitnesses, *p);
			break;
		case 8:
			for (int i = 0; i < numofparticles; i++)
				{
				Xmax[i] = 500;
				Xmin[i] = -500;
				}
			F8(X, fitnesses, *p);
			break;
		case 9:
			F9(X, fitnesses, *p);
			break;
		case 10:
			F10(X, fitnesses, *p);
			break;
		case 11:
			F11(X, fitnesses, *p);
			break;
		case 12:
			F12(X, fitnesses, *p);
			break;
		case 13:
			F13(X, fitnesses, *p);
			break;
		case 14:
			F14(X, fitnesses, *p);
			break;
		case 15:
			F15(X, fitnesses, *p);
			break;
		default:
			;
		}

    minfit = *min_element(fitnesses, fitnesses + numofparticles);
    minfitidx = min_element(fitnesses, fitnesses + numofparticles) - fitnesses;
    
		// STEP 6. set gbest : The best position of the swarm
		if(minfit < *gbestfit){
            *gbestfit = minfit;
            memcpy(gbest, X[minfitidx], sizeof(float) * numofdims);
            cout << t << " times fitness is " << minfit << endl;
        }
        
        bests[t] = *gbestfit;

    }
}


int main(){
	
    time_t t;
    srand((unsigned) time(&t));

    float xmin[numofparticles], xmax[numofparticles]; // value declaration
    float initpop[numofparticles][numofdims];
    float bests[1000];
    float gbestfit, point;
    float gbest[numofdims];
    int select;
    
    for(int i = 0; i < 10; i++){ // set lower and upper value
        xmax[i] = 100;
        xmin[i] = -100;
    	}
    	
    for(int i = 0; i < 20; i++) //set all particles to random values
        for(int j = 0; j < 10; j++){ 
           	initpop[i][j] = rand() % (100 + 100 + 1) - 100;
        }
    
    
    while(1){ // Select the Questions.
	  cout << endl << "         Final Project - Park Hoijai (18309002)" << endl;
	  cout << endl << "=======================================================" << endl << endl;
	  cout << "01. Question 1 - Sphere function Equation" << endl;
	  cout << "02. Question 2 - Schwefel's 2.22 function Equation" << endl;
	  cout << "03. Question 3 - Schwefel's 1.20 function Equation" << endl;
	  cout << "04. Question 4 - Schwefel's 2.21 function Equation" << endl;
	  cout << "05. Question 5 - Rosenbrock function Equation" << endl;
	  cout << "06. Question 6 - Step function Equation" << endl;
	  cout << "07. Question 7 - Quartic noise function Equation" << endl;
	  cout << "08. Question 8 - Schwefel function Equation" << endl;
	  cout << "09. Question 9 - Rastrigin function Equation" << endl;
	  cout << "10. Question 10 - Ackley function Equation" << endl;
	  cout << "11. Question 11 - Griewank function Equation" << endl;
	  cout << "12. Question 12" << endl;
	  cout << "13. Question 13" << endl;
	  cout << "14. Question 14" << endl;
	  cout << "15. Question 15" << endl << endl;	  
	  cout << "=======================================================" << endl;
	  cout << "Command : ";
	  
	  cin >> select;
	  cout << endl;
	  	if(select>15 || select<=0){
			cout<<"Program is terminated. Pleas insert number 1~15. \n";
			exit;
		}else{
			PSO(1000, 2, 2, xmin, xmax, initpop, bests, &gbestfit, gbest, &point, select);
			cout << endl << "Result : " << gbestfit << endl;
	  		cout << endl << "================== Question." << select << " Finish ==================" << endl;
		}
	}	
	return 0;
}
