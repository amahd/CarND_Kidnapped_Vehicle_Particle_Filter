/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

#define MAX_PART 20
default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	//  Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


	num_particles = MAX_PART; //change macro to control particle number


	Particle temp;




	// defining normal dists for coordinates (x,y,theta), using current location as mean
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);

	// Now initializing all particles with random position
	for (int h =0; h < num_particles; ++h)  {

		temp.id = h;                       	//ID
		temp.x = dist_x(gen);			  	// Noisy X
		temp.y = dist_y(gen);				// Noisy y
		temp.theta = dist_psi(gen);			// Noisy angle
		temp.weight = 1;					// weight of particle
		particles.push_back(temp);
	}

	// Initialize also all weights to 1 for variable weights
	for (int h =0; h < num_particles; ++h)
		weights.push_back(1);


	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// : Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/




	double x,y,th;   // new position coordinates

	for (auto it = particles.begin(); it!=particles.end(); ++it) {

		// Get new x,y,theta from previous location and speed

		if (std::abs(yaw_rate) > 0.0001){
			x =  predx(it->x, it->theta,velocity, yaw_rate,delta_t);
			y =  predy(it->y, it->theta,velocity, yaw_rate,delta_t);
			th = predth(it->theta,yaw_rate,delta_t );
		 }
		else {
			th = it->theta;
			x =   it->x + velocity*cos(it->theta)*delta_t;
			y =   it->y + velocity*sin(it->theta)*delta_t;


		}



		//  temper the above values with Normal noise
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_psi(th, std_pos[2]);

		// Store the predicted noisy location in the particle
		it->x = dist_x(gen);
		it->y = dist_y(gen);
		it->theta = dist_psi(gen);

	}



}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	//  Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// : Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution



	//Go through each particle and scan for landmkarks and observations


	for (auto particle = particles.begin(); particle!=particles.end(); ++particle) { // All particles
		Map range_map;


		for (int j = 0; j < map_landmarks.landmark_list.size(); ++j){  // loop for very landmarks to be in sensor range
			// if any land mark is in sensor range of car, save it
			if (dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particle->x,particle->y ) <= sensor_range)
				range_map.landmark_list.push_back(map_landmarks.landmark_list[j]);

		    }


		// Convert all observations per particle to global coordinates
		double xcor,ycor;	//x and y coordinates for car->global conversion
		std::vector<LandmarkObs> global_obs;  //all obersations in global coordinates
		LandmarkObs temp;

		for (auto ob = observations.begin(); ob!=observations.end(); ++ob) {
			xcor = particle->x + ob->x*cos(particle->theta)- ob->y*sin(particle->theta);    // get x-coordinate of landmark in map domain
			ycor = particle->y + ob->y*cos(particle->theta)+ ob->x*sin(particle->theta);	// get y-coordinate of landmark in map domain
			temp.id = -1;  // invalid ID
			temp.x = xcor;
			temp.y = ycor;
			global_obs.push_back(temp);
		    }

		// Compare all observations with global map and find nearest neighbour

		for (auto ob = global_obs.begin(); ob!=global_obs.end(); ++ob){  //check for all observations
			double min_dist = numeric_limits<double>::infinity(); // assign a high value and look for smaller every time
			double actual_dist;
				for (int j = 0; j < range_map.landmark_list.size(); ++j){ //loop through all landmakrs in map

					actual_dist = dist(ob->x, ob->y, range_map.landmark_list[j].x_f, range_map.landmark_list[j].y_f );
					if (actual_dist < min_dist){
						min_dist = actual_dist;
						ob->id = range_map.landmark_list[j].id_i;
					   }

				   }
			  }


		particle->weight = 1.;
		double std_x = std_landmark[0];
	    double std_y = std_landmark[1];
		double std_x_2 = 2 * std_x * std_x;
		double std_y_2 = 2 * std_y * std_y;
		double denom = 2 * M_PI * std_x * std_y;


		//Valid Observations are used to calculate weights

		for (unsigned j = 0; j < global_obs.size(); j++) {
			LandmarkObs obs = global_obs[j];

			//use intermediate variable for product to ease debug
			double ans =  exp(-pow(obs.x-map_landmarks.landmark_list[obs.id-1].x_f, 2)/std_x_2 -pow(obs.y - map_landmarks.landmark_list[obs.id-1].y_f, 2)/std_y_2)/denom;
			particle->weight *=ans;
		}

	}

	for (int h =0; h < num_particles; ++h)
		weights[h] = particles[h].weight;


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution





		vector<Particle> final_particles;    // final selected particles
		default_random_engine gen;
		discrete_distribution<> distribution (weights.begin(), weights.end());

		for (int i = 0; i < num_particles; ++i)	{
			int idx = (int) distribution(gen);
			final_particles.push_back(particles[idx]);  //move particles with highest probabilities to final
		}


		//overwrite particles list
		particles = final_particles;


		//reset weights to 1
		for (int h =0; h < num_particles; ++h)
			weights[h] = 1;




}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
