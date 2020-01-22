//
//  main.cpp
//  Assignment2
//
//  Created by Sultanmurat on 10/4/19.
//  Copyright © 2019 Sultanmurat. All rights reserved.
//
#include<iostream>
#include<cmath>
#include<random>
#include<fstream>
#include<cmath>
#include<array>
#include"/Users/zhazira/Desktop/Evolution/Assignment1/Assignment1/gnuplot_i.hpp"
using namespace std;

std::random_device rd; // obtain a random number from hardware
std::mt19937 eng(rd()); // seed the generator
std::uniform_real_distribution<> distr01(0.0, 1.0);
std::uniform_real_distribution<> distr02(0.0, 1.0);

double f_x(double x);
void create_points(int n, int low, int high);
void plot_points();
void wait_for_key();
void plot_equation();
double x, y;
struct Point{
    double x_coord;
    double y_coord;
};
double crossover_probability = 1;
double errors[10000000] = {0};
double best_fitnesses[10000000] = {0};
vector<Point> points;
vector<double> errorsHC;
// * 2, / 3, + 4, - 5, sin 6, cos 7, constant 0, x 1
const int number_population = 100;
void read_points();
vector<double> mutateConstant(vector<double> child1);
vector<double> mutateAddTree(vector<double> child1);
bool compareIndividuals(vector<double> ind1, vector<double> ind2);
vector<vector<double>> crossover(vector<double> parent1, vector<double> parent2);
vector<vector<double>> crossover2(vector<double> parent1, vector<double> parent2);
vector<double> create_tree_grow(int max_depth);
vector<double> create_tree_full(int max_depth);
double calculate_function(double x, vector<double> heap);
int find_depth(vector<double> heap){
    int max = 1;
    for (int i = 1; i < heap.size(); i++){
        if (heap[i] != 99) max = i;
    }
    int depth = (int) log2(max) + 1;
    return depth;
}
double calculate_fitness2(vector<double> heap);
double calculate_fitness(vector<double> heap){
    size_t n_points = points.size();
     double res = 0;
     for (int i = 0; i < n_points; i++){
         double x = points[i].x_coord;
         double y_exp = calculate_function(x, heap);
         double y_real = points[i].y_coord;
         res += abs(y_exp-y_real);
     }
     res = res/n_points;
     return res;
}
int complexity(vector<double> heap);
ofstream dataCutted;
vector<double> mutateAddConstant(vector<double> heap);
vector<double> mutateMultiplyConstant(vector<double> heap);
vector<double> snipping(vector<double> heap);
vector<double> pruning(vector<double> heap);
void drawTree(vector<double> heap);
int main(int argc, const char * argv[]) {
    //dataCutted.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/dataCutted.txt");
    vector<vector<double>> population(number_population);
    create_points(1000, 0, 10);
    read_points();
//    ofstream learningCurveHC2;
//    learningCurveHC2.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveHC2.txt");
//    for (int i = 0 ; i < points.size(); i ++){
//        if (points[i].x_coord <= 50000 ){
//            if (points[i].y_coord > 1) points[i].y_coord = 1;
//            if (errorsHC[i] > 1) errorsHC[i] = 1;
//            learningCurveHC2 << points[i].x_coord/100 << " " << points[i].y_coord << " " << errorsHC[i] << endl;
//        }
//    }
   // wait_for_key();
//
//    for (int i = 0; i < 250; i++){
//        int current = (int)(distr01(eng)*1000);
//        dataCutted << points[current].x_coord << " " << points[current].y_coord << endl;
//    }
   // wait_for_key();
    size_t n_points = points.size();
//    for (int i = 0; i < 20; i++){
//        int max = 6;
//        int r1 =(int) (distr01(eng)*(pow(2,max-1)-1))+1;
//        int r2 = (int) (distr01(eng)*(pow(2, (int) log2(r1)))) + pow(2, (int) log2(r1));
//        if (r2 == r1 && r2 != 1){
//            r2 = (pow(2, (int)log2(r1)+1)-1) - (r1 - pow(2, (int)log2(r1)));
//        }
//              cout << "r1: " << r1 << "; r2: " << r2 << endl;
//    }
   // code to check snipping
//    vector<double> myTree(256);
//    for (int i = 0; i < myTree.size(); i++){
//        myTree[i] = 99;
//    }
//    myTree[1] = 106;
//    myTree[2] = 200.523;
//    myTree[3] = 103;
//    myTree[4] = 99;
//    myTree[6] = 101;
//    myTree[7] = 106;
//    myTree[8] = 99;
//    myTree[14] = 101;

//    vector<double> mutated = snipping(myTree);
//    cout << "My tree: ";
//    for (int i = 0; i < 16; i++){
//        cout << mutated[i] << endl;
//    }
//    mutated = pruning(mutated);
//    cout << endl << "Mutated tree: ";
//       for (int i = 0; i < 256; i++){
//           cout << mutated[i] << endl;
//       }
//
//    drawTree(mutated);
//
//    wait_for_key();
    
     //create initial population
    int twos = 0;
    int j = 2;
    for (int i = 0; i < number_population; i++){
        //cout << j << endl;
        //population[i] = create_tree_grow(j);
        //if (distr01(eng) > 0.5){
            population[i] = create_tree_full(j);
        //} else{
        //   population[i] = create_tree_grow(j);
        //}
                //cout << "Fitness1: " << calculate_fitness(population[i]) << endl;
        //cout << "Fitness2: " << calculate_fitness2(population[i]) << endl;
        //cout << "Depth created: " << j << "; Depth found: " << find_depth(population[i]) << endl;
        if (j==2) twos++;
        j++;
        if (j>3) {
            if (twos > 3){
                j = 3;
            } else {
                j=2;
            }
        }
    }
    j = 0;
    int number_generations = 200;
    //wait_for_key();
    
    
    // code to check crossover
//    population[0] = create_tree_full(8);
//    population[1] = create_tree_full(8);
//    for (int i = 0;i <16; i++){
//        cout << population[0][i] << endl;
//
//    }
//    cout << endl;
//
//    for (int i = 0;i < 16; i++){
//        cout << population[1][i] << endl;
//    }
//    cout << endl;
//
//    vector<vector<double>> children = crossover2(population[0], population[1]);
//    vector<double> child1 = children[0];
//    vector<double> child2 = children[1];
//    for (int i = 0;i < 16; i++){
//        cout << child1[i] << endl;
//    }
//    cout << endl;
//
//    for (int i = 0;i < 16; i++){
//        cout << child2[i] << endl;
//    }
//
//    cout << endl;
//    wait_for_key();
 //code to check mutation
    //population[0] = create_tree_full(3);
//    vector<double> randTree(256);
//    for (int i = 0; i < 256; i++){
//        randTree[i] = 99;
//    }
//    randTree[1] = 105;
//    randTree[2] = 101;
//    randTree[3] = 101;
//    cout << endl << endl << "Original tree: " << endl;
//       for (int i = 1; i < 16; i++){
//           cout << randTree[i] << endl;
//       }
//       vector<double> mutated = snipping(randTree);
//       cout << endl << endl << "Mutated tree: "<< endl;
//       for (int i = 1; i < 16; i++){
//           cout << mutated[i] << endl;
//       }
//    wait_for_key();
    // Random search
//    vector<vector<double>> min_trees;
//    int index_min = 0;
//    double min = 10;
//    int k = 2;
//    ofstream randomCurve;
//    randomCurve.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/randomCurve.txt");
//
//    while (j <= number_generations){
//        k = (int) (distr01(eng)*6) + 3;
//        vector<double> current_tree = create_tree_full(3);
////        if (k>8) {
////            if (twos < 6) {
////                k = 2;
////                twos+=1;
////            }
////            else k = 3;
////        }
//        double current_fitness = calculate_fitness2(current_tree);
//        if (current_fitness < min){
//            min = current_fitness;
//            min_trees.push_back(current_tree);
//            index_min = j;
//        }
//        best_fitnesses[j] = min;
//        j++;
//        randomCurve << j << " " << min << endl;
//        cout << j << " " << min << endl;
//    }
//    ofstream randomEquation;
//    randomEquation.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/randomEquation.txt");
//    vector<double> best_tree = min_trees[min_trees.size()-1];
//
//    for (int i = 0; i < n_points; i++){
//        randomEquation << points[i].x_coord << " "<<calculate_function(points[i].x_coord, best_tree) << endl;
//    }
//    Gnuplot gpCurve("lines");
//    gpCurve.set_title("Learning Curve for Random Search");
//    gpCurve.set_xlabel("Evaluations");
//    gpCurve.set_ylabel("MAE Error");
//    gpCurve.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/randomCurve.txt", 1, 2, "Learning curve of Random Search");
//    gpCurve.showonscreen();
//
//    Gnuplot gp("points");
//    gp.set_xlabel("x");
//    gp.set_ylabel("y");
//    gp.set_title("The best equation found");
//    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/data.txt", 1, 2, "Data points");
//    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/randomEquation.txt", 1, 2, "Equation");
//    gp.showonscreen();
//    for (int i = 0; i< best_tree.size(); i++){
//        cout << best_tree[i] << endl;
//    }
//    wait_for_key();
    
    // hill climber
//    ofstream learningCurveHC;
//    learningCurveHC.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveHC.txt");
//    vector< vector<double> > current_trees(number_population);
//    for (int i = 0; i < number_population; i++){
//        int r = (int) (distr01(eng)*4) + 3;
//        //cout << "r" << r << endl;
//        vector<double> current_tree = create_tree_full(3);
//        current_trees[i] = current_tree;
//    }
//    while (j < number_generations  ){
//        double current_fitnesses[number_population];
//        double total_fitnesses = 0;
//        for (int i = 0; i < number_population; i++){
//            if (j+1%99 == 0) current_trees[i] = snipping(current_trees[i]);
//            if (j+1%100==0) current_trees[i] = pruning(current_trees[i]);
//
//            double current_fitness = calculate_fitness2(current_trees[i]);
//            total_fitnesses += current_fitness;
//            current_fitnesses[i] = current_fitness;
//            vector<double> mutated;
//            double randomNumber = distr01(eng);
//            if (randomNumber > 0.75 ){
//                mutated = mutateAddConstant(current_trees[i]);
//            } else if (randomNumber > 0.5 ){
//                mutated = mutateMultiplyConstant(current_trees[i]);
//            } else if (randomNumber > 0.25 ){
//                mutated = mutateConstant(current_trees[i]);
//            } else {
//                mutated = mutateAddTree(current_trees[i]);
//            }
//            double mutated_fitness = calculate_fitness2(mutated);
//            if (mutated_fitness < current_fitness) current_trees[i] = mutated;
//        }
//
//
//        double mean = total_fitnesses/number_population;
//        double error_deviation = 0;
//        for (int i = 0; i < number_population; i++){
//            error_deviation = error_deviation + sqrt((current_fitnesses[i]-mean)*(current_fitnesses[i]-mean));
//        }
//        error_deviation = error_deviation/number_population;
//        cout << j << " " << mean << " " << error_deviation << endl;
//        learningCurveHC << j << " " << mean << " " << error_deviation << endl;
//        j++;
//    }
//    vector<double> hc_best_tree = current_trees[0];
//    drawTree(hc_best_tree);
//    ofstream HCEquation;
//    HCEquation.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/HCEquation.txt");
//    for (int i = 0; i < n_points; i++){
//        HCEquation << points[i].x_coord << " "<< calculate_function(points[i].x_coord, hc_best_tree) << endl;
//    }
//    Gnuplot hc("points");
//    hc.set_xlabel("x");
//    hc.set_ylabel("y");
//    hc.set_title("The best equation found using HC");
//    hc.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/data.txt", 1, 2, "Data points");
//    hc.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/HCEquation.txt", 1, 2, "Equation");
//    hc.showonscreen();
//
//    Gnuplot HCCurve("lines");
//    HCCurve.set_title("Learning Curve for HC Search");
//    HCCurve.set_xlabel("Evaluations");
//    HCCurve.set_ylabel("MAE Error");
//HCCurve.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveHC.txt", 1, 2, 3, "Learning curve for HC Search");
//      HCCurve.showonscreen();
//
//
//
//    wait_for_key();

    // Genetic Programming

    double min;
    ofstream learningCurveGP2;
    ofstream diversity;
    ofstream dot;
    ofstream complexityFile;
    complexityFile.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/complexity.txt");
    dot.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/dot.txt");
    diversity.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/diversity.txt");
    learningCurveGP2.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveGP2.txt");
    while (j < number_generations){
        vector<vector<double>> next_population;
        // Get fitness values
        array<double, 2> fitness[number_population];
        double total_fitness = 0;
        double max_fitness = 0;
        double first_fitness = 10000;
        double second_fitness = 10000;
        int first_index = 0;
        int second_index = 0;
        if (j+1%4==0){
            for (int i = 0; i < number_population; i++){
                population[i] = snipping(population[i]);
            }
        }
        if (j+1%5==0){
            for (int i = 0; i < number_population; i++){
                population[i] = pruning(population[i]);
            }
        }
        for (int i = 0; i < number_population; i++){
            double current_fitness = calculate_fitness2(population[i]);
            int complexityNumber = complexity(population[i]);
            if (current_fitness < first_fitness){
                second_fitness = first_fitness;
                second_index = first_index;
                first_fitness = current_fitness;
                first_index = i;
            }
            if (current_fitness > first_fitness && current_fitness < second_fitness){
                second_fitness = current_fitness;
                second_index = i;
            }
            if (current_fitness > max_fitness){
                max_fitness = current_fitness;
            }
            fitness[i][0] = current_fitness;
            fitness[i][1] = i;
            if (current_fitness < 1){
                dot << j << " " << current_fitness << endl;
                complexityFile << current_fitness << " " << complexityNumber << endl;
            }
                //cout << "i: " << i << "  fitness: "  << current_fitness << endl;
            total_fitness += fitness[i][0];
        }

        // find error
        double mean = total_fitness/number_population;
        double error_deviation = 0;
        for (int i = 0; i < number_population; i++){
            error_deviation = error_deviation + sqrt((fitness[i][0]-mean)*(fitness[i][0]-mean));
        }
        error_deviation = error_deviation/sqrt(number_population);
        error_deviation = error_deviation/number_population;
        errors[j] = error_deviation;
        if (error_deviation > 0.5) error_deviation = 0.5;

        if (j+1%17==0){
            for (int i = 0; i < number_population; i++){
                population[i] = snipping(population[i]);
            }
        }
        if (j+1%18==0){
            for (int i = 0; i < number_population; i++){
                population[i] = pruning(population[i]);
            }
        }

        // Producing new children
        next_population.push_back(population[second_index]);
        next_population.push_back(population[first_index]);
        for (int m =0; m < number_population-2; m++){

            //Select using tournament selection
            double k = number_population/30;  // set selection pressure
            int r1 = -1;//(int) (distr01(eng)*number_population);
            int r2 = -1; //(int) (distr01(eng)*number_population);
            for (int i =0; i < k; i++){
                int ind = (int) (distr01(eng)*number_population);
                if (r1 == -1) {
                    r1 = ind;
                }else {
                    double current = fitness[ind][0];
                    double best = fitness[r1][0];
                    //cout << ind << endl;
                    if (current < best){
                        r1 = ind;
                    }
                }
            }
            for (int i =0; i < k; i++){
                int ind = (int) (distr01(eng)*number_population);
                if (r2 == -1) {
                    r2 = ind;
                }else {
                    double current = fitness[ind][0];
                    double best = fitness[r2][0];
                    //cout << ind << endl;
                    if (current < best){
                        r2 = ind;
                    }
                }
            }
            //cout << "r1: " << r1 << "   r2: " << r2 << endl;


            vector<double> parent1 = population[r1];
            vector<double> parent2 = population[r2];

            // pushing children values into new population
            vector<vector<double>> children;
            children = crossover2(parent1, parent2);
            vector<double> child1;
            vector<double> child2;
            child1 = children[0];
            child2 = children[1];
            child1 = mutateAddTree(child1);
            child2 = mutateAddTree(child2);
            child1 = mutateConstant(child1);
            child2 = mutateConstant(child2);


            child1 = mutateAddConstant(child1);
            child2 = mutateAddConstant(child2);

            child1 = mutateMultiplyConstant(child1);
            child2 = mutateMultiplyConstant(child2);






            double fCh1 = calculate_fitness2(child1);
            double fCh2 = calculate_fitness2(child2);
            double dp1c1 = abs(fitness[r1][0] - fCh1);
            double dp2c2 = abs(fitness[r2][0] - fCh2);
            double dp1c2 = abs(fitness[r1][0] - fCh2);
            double dp2c1 = abs(fitness[r2][0] - fCh1);
            if (fCh1 < fCh2){
                next_population.push_back(child1);
            } else {
                next_population.push_back(child2);
            }
//            if (dp1c1 + dp2c2 < dp1c2+dp2c1){
//                // compare c1 to p1
//                if (fCh1 < fitness[r1][0]){
//                    next_population.push_back(child1);
//                } else {
//                    next_population.push_back(parent1);
//                }
//                if (fCh2 < fitness[r2][0]){
//                    next_population.push_back(child2);
//                } else {
//                    next_population.push_back(parent2);
//                }
//            } else {
//                if (fCh1 < fitness[r2][0]){
//                    next_population.push_back(child1);
//                } else {
//                    next_population.push_back(parent2);
//                }
//                if (fCh2 < fitness[r1][0]){
//                    next_population.push_back(child2);
//                } else {
//                    next_population.push_back(parent1);
//                }
//            }
        }

        population = next_population;
        j++;
        double diversityNumber = (max_fitness - first_fitness)/number_population;
        //if (diversityNumber > 5) diversityNumber = 5;

        cout << "Gen: " << j << "   Min: " << setprecision(5) << first_fitness << "   Diff: " << setprecision(5) << diversityNumber
                <<endl;

        // Plots
        learningCurveGP2 << j << " " << first_fitness  << " " << error_deviation << endl;
        diversity << j << " " << diversityNumber << endl;

        if (diversityNumber < 0.0000000001 ) break;
    }

    min = calculate_fitness2(population[0]);
    int index_min = 0;
    for (int i = 0; i < number_population; i++){
        double current_fitness = calculate_fitness2(population[i]);
        if (current_fitness < min ){
            min = current_fitness;
            index_min = i;
        }
    }
    vector<double> best_tree = population[index_min];
    best_tree = snipping(best_tree);
    //best_tree = pruning(best_tree);


    cout << "The best equation found: " << endl;
    for (int i = 1; i < 256; i++){
        cout << best_tree[i] << endl;
    }

    ofstream GPEquation;
    GPEquation.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/GPEquation.txt");


    drawTree(best_tree);
    for (int i = 0; i < n_points; i++){
        GPEquation << points[i].x_coord << " "<< calculate_function(points[i].x_coord, best_tree) << endl;
    }

    Gnuplot gpCurve("lines");
    gpCurve.set_title("Learning Curves");
    gpCurve.set_xlabel("Evaluations (x 100)");
    gpCurve.set_ylabel("MAE Error");
    //gpCurve.set_ylogscale();
   gpCurve.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveGP2.txt", 1, 2,3, "Learning curve for GP Search v2");
    //gpCurve.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveHC2.txt", 1, 2,3, "Learning curve for HC Search");
    //gpCurve.plotfile_xy_err("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/learningCurveGP.txt", 1, 2,3, "Learning curve for GP Search v1");
    
    //gpCurve.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/randomCurve2.txt", 1, 2, "Learning curve for Random Search");
      gpCurve.showonscreen();
//    wait_for_key();

    Gnuplot gp("points");
    gp.set_xlabel("x");
    gp.set_ylabel("y");
    gp.set_title("The best equation found");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/points.txt", 1, 2, "Data points");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/GPEquation.txt", 1, 2, "Equation");
    gp.showonscreen();

    Gnuplot gpDot("dots");
    gpDot.set_title("Dot for GP");
    gpDot.set_xlabel("Generation");
    gpDot.set_ylabel("Fitness");

    gpDot.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/dot.txt", 1, 2);
    gpDot.showonscreen();

    Gnuplot gpDiversity("lines");
    gpDiversity.set_title("Diversity for GP");
    gpDiversity.set_xlabel("Generation");
    gpDiversity.set_ylabel("Fitness");
    gpDiversity.set_ylogscale();
    gpDiversity.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/diversity.txt", 1, 2, "Diversity for GP");
    gpDiversity.showonscreen();

    Gnuplot gpComplexity("points");
    gpComplexity.set_title("Complexity vs Accuracy");
    gpComplexity.set_xlabel("Fitness value");
    gpComplexity.set_ylabel("Complexity");

    gpComplexity.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/complexity.txt", 1, 2);
    gpComplexity.showonscreen();

      wait_for_key();
    
    return 0;
}

void drawTree(vector<double> heap){
    ofstream outputTree;
    outputTree.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/outputTree.txt");
    int depth = find_depth(heap);
    int max_cols = pow(2,depth);
    for (int i = depth; i >=1; i--){
        int n_nodes = pow(2, i-1);
        int area_for_node = max_cols/n_nodes;
        for (int j = 1; j <= n_nodes; j++){
            for (int k = 1; k <= area_for_node; k++){
                if (k!=area_for_node/2) outputTree << " ";
                else if (k==area_for_node/2){
                    switch ((int)heap[pow(2,i-1)+j-1])
                    {
                        
                        case 99: outputTree << " ";
                            break;
                        case 101: outputTree << "x";
                            break;
                        case 102: outputTree << "*";
                            break;
                        case 103: outputTree << "/";
                            break;
                        case 104: outputTree << "+";
                            break;
                        case 105: outputTree << "-";
                            break;
                        case 106: outputTree << "s";
                            break;
                        case 107: outputTree << "c";
                            break;
                        default:
                            if (heap[pow(2,i-1)+j-1]-200 < 0) {
                                outputTree << std::fixed;
                                outputTree << setprecision(1) << heap[pow(2,i-1)+j-1]-200 << " " ;
                            } else{
                                outputTree << std::fixed;
                                outputTree << setprecision(1) << heap[pow(2,i-1)+j-1]-200;
                            }
                                break;
                    }
                }
            }
        }
        outputTree << endl;
    }
}

int complexity(vector<double> heap){
    int complexity = 0;
    size_t n = heap.size();
    for (int i = 0; i < n; i++){
        if (heap[i] ==  101) complexity+=2;
        else if (heap[i] >101 && heap[i] < 106) complexity+=3;
        else if (heap[i] == 106 || heap[i] == 107) complexity +=4;
        else complexity+=1;
    }
    return complexity;
}

double calculate_fitness2(vector<double> heap){
    size_t n_points = points.size();
    size_t n = heap.size();
    double res = 0;
    double res_heap[n_points][n];
    for (int i = n-1; i >=1; i--){
        for (int j = 0; j <n_points; j++){
            res_heap[j][i] = heap[i];
        }
    }
    for (int i = n-1; i >= 1; i--){
        for (int j = 0; j < n_points; j++){
            double x = points[j].x_coord;
            switch ((int)heap[i])
            {
                case 99:
                    break;
                case 101: res_heap[j][i] = x;
                    break;
                case 102: res_heap[j][i] = res_heap[j][i*2]*res_heap[j][i*2+1];
                    break;
                case 103: if (res_heap[j][i*2+1] != 0) res_heap[j][i] = res_heap[j][i*2]/res_heap[j][i*2+1];
                else res_heap[j][i] = 1;
                    break;
                case 104: res_heap[j][i] = res_heap[j][i*2]+res_heap[j][i*2+1];
                    break;
                case 105: res_heap[j][i] = res_heap[j][i*2]-res_heap[j][i*2+1];
                    break;
                case 106: res_heap[j][i] = sin(res_heap[j][i*2]);
                    break;
                case 107: res_heap[j][i] = cos(res_heap[j][i*2]);
                    break;
                default: res_heap[j][i] = res_heap[j][i]-200.0;
                    //cout << "const: " << res[i] << endl;
                    break;
            }
        }
    }
    for (int j = 0; j < n_points; j++){
        double y_exp = points[j].y_coord;
        double y_real = res_heap[j][1];
        res += abs(y_exp-y_real);
    }
    res = res/n_points;
    return res;
}

vector<double> snipping(vector<double> heap){
    size_t n_points = points.size();
    size_t n = heap.size();
    vector<double> mutated(n);
    for (int i = 0; i < n; i++){
        mutated[i] = heap[i];
    }
    double res_heap[n_points][n];
    for (int i = n-1; i >=1; i--){
        for (int j = 0; j <n_points; j++){
            res_heap[j][i] = heap[i];
        }
    }
    
    for (int i = n-1; i >= 1; i--){
        double min = 10;
        double max = 0;
        for (int j = 0; j < n_points; j++){
            double x = points[j].x_coord;
            switch ((int)heap[i])
            {
                case 99:
                    break;
                case 101: res_heap[j][i] = x;
                    break;
                case 102: res_heap[j][i] = res_heap[j][i*2]*res_heap[j][i*2+1];
                    break;
                case 103: if (res_heap[j][i*2+1] != 0) res_heap[j][i] = res_heap[j][i*2]/res_heap[j][i*2+1];
                else res_heap[j][i] = 1;
                    break;
                case 104: res_heap[j][i] = res_heap[j][i*2]+res_heap[j][i*2+1];
                    break;
                case 105: res_heap[j][i] = res_heap[j][i*2]-res_heap[j][i*2+1];
                    break;
                case 106: res_heap[j][i] = sin(res_heap[j][i*2]);
                    break;
                case 107: res_heap[j][i] = cos(res_heap[j][i*2]);
                    break;
                default: res_heap[j][i] = res_heap[j][i]-200.0;
                    //cout << "const: " << res[i] << endl;
                    break;
            }
            if (heap[i] > 101 && heap[i] < 108){
                if (j==0) {
                    min = res_heap[j][i];
                    max = res_heap[j][i];
                }
                if (res_heap[j][i] < min){
                    min = res_heap[j][i];
                }
                if (res_heap[j][i] > max){
                    max = res_heap[j][i];
                }
            }
        }
        if (abs(max - min) < 0.4) {
            //cout << "Max " << max << "  Min " << min << endl;
            //cout << "i: " << i << endl;
            int power = 0;
            for (int m = i; m < 256; m=m*2){
                for (int l = 0; l<pow(2, power); l++){
                    if (m==i){
                        heap[i] = 200+(max+min)/2;
                        mutated[i] = 200+ (max+min)/2;
                    } else{
                        mutated[m+l] = 99;
                        heap[m+l] =99;
                    }
                }
                power++;
            }
        }
//        } else if (abs(abs(max)-abs(min))<0.001  && abs(max) < 0.9){
//            int power = 0;
//            for (int m = i; m < 256; m=m*2){
//                for (int l = 0; l<pow(2, power); l++){
//                    if (m==i){
//                        heap[i] = 106;
//                        mutated[i] = 106;
//                    } else if (m==i*2 && l==0){
//                        heap[m+l] = 101;
//                        mutated[m+l] = 101;
//                    }else{
//                        mutated[m+l] = 99;
//                        heap[m+l] =99;
//                    }
//                }
//                power++;
//            }
//        }
    }
        
    return mutated;
}


vector<double> pruning(vector<double> heap){
    size_t n_points = points.size();
    size_t n = heap.size();
    vector<double> mutated(n);
    for (int i = 0; i < n; i++){
        mutated[i] = heap[i];
    }
    double offset = 0.1;
    for (int i = 1; i <n; i++){
        if (mutated[i] == 102){
            if (abs(mutated[i*2]-200) < 1+offset && abs(mutated[i*2]-200) > 1-offset){
                
                mutated[i*2] = 99;
                int power = 0;
                int k = i*2+1;
                for (int m = i; m < 256; m=m*2){
                    for (int l = 0; l<pow(2, power); l++){
                        mutated[m+l] = heap[k+l];
                    }
                    k = k*2;
                    power++;
                }
            } else if (abs(mutated[i*2+1]-200) < 1+offset && abs(mutated[i*2+1]-200) > 1-offset){
               
                mutated[i*2+1] = 99;
                int power = 0;
                int k = i*2;
                for (int m = i; m < 256; m=m*2){
                    for (int l = 0; l<pow(2, power); l++){
                        mutated[m+l] = heap[k+l];
                    }
                    k = k*2;
                    power++;
                }
            }
        } else if (mutated[i] == 104){
            if (abs(mutated[i*2]-200) < offset){
                mutated[i*2] = 99;
                int power = 0;
                int k = i*2+1;
                for (int m = i; m < 256; m=m*2){
                    for (int l = 0; l<pow(2, power); l++){
                        mutated[m+l] = heap[k+l];
                    }
                    k = k*2;
                    power++;
                }
            } else if (abs(mutated[i*2+1]-200) < offset ){
                mutated[i*2+1] = 99;
                int power = 0;
                int k = i*2;
                for (int m = i; m < 256; m=m*2){
                    for (int l = 0; l<pow(2, power); l++){
                        mutated[m+l] = heap[k+l];
                    }
                    k = k*2;
                    power++;
                }
            }
        } else if (mutated[i]== 103){
            if (abs(mutated[i*2+1]-200) <1+offset && abs(mutated[i*2+1]-200) > 1-offset){
                
                mutated[i*2+1] = 99;
                int power = 0;
                int k = i*2;
                for (int m = i; m < 256; m=m*2){
                    for (int l = 0; l<pow(2, power); l++){
                        mutated[m+l] = heap[k+l];
                    }
                    k = k*2;
                    power++;
                }
            }
        } else if (mutated[i] == 105){
            if (abs(mutated[i*2+1]-200) < 0.1 ){
                mutated[i*2+1] = 99;
                int power = 0;
                int k = i*2;
                for (int m = i; m < 256; m=m*2){
                    for (int l = 0; l<pow(2, power); l++){
                        mutated[m+l] = heap[k+l];
                    }
                    k = k*2;
                    power++;
                }
            }
        }
    }
    return mutated;
}

vector<double> mutateAddConstant(vector<double> heap){
    if (find_depth(heap) < 8){
        double r = distr01(eng)*0.3;
        vector<double> res(256);
        for (int i= 0; i < res.size(); i++){
            res[i] = 99;
        }
        res[1] = 104;
        if (distr01(eng) > 0.5){
            res[2] = 200+r;
        } else {
            res[2] = 200-r;
        }
        int power = 0;
        int k = 1;
        for (int i = 3; i < 256; i=i*2){
            for (int l = 0; l<pow(2, power); l++){
                res[i+l] = heap[k+l];
            }
            k = k*2;
            power++;
        }
        return res;
    } else {
        //int r = (int) (distr01(eng)*128) + 128;
        //heap[r] = 101;
        return heap;
    }
}

vector<double> mutateMultiplyConstant(vector<double> heap){
    if (find_depth(heap) < 8){
        double r = distr01(eng)*0.3+1;
        vector<double> res(256);
        for (int i = 0; i < res.size(); i++){
            res[i] = 99;
        }
        res[1] = 102;
        if (distr01(eng) > 0.5){
            res[2] = 200+r;
        } else {
            res[2] = 200-r;
        }
        int power = 0;
        int k = 1;
        for (int i = 3; i < 256; i=i*2){
            for (int l = 0; l<pow(2, power); l++){
                res[i+l] = heap[k+l];
            }
            k = k*2;
            power++;
        }
        return res;
    } else {
        //int r = (int) (distr01(eng)*128) + 128;
        //heap[r] = 101;
        return heap;
    }
}

vector<double> mutateConstant(vector<double> child1){
    size_t n_nodes = child1.size();
    vector<double> temp_child = child1;
    for (int i = 0; i < n_nodes; i++){
        if (child1[i] <= 210 && child1[i]>= 190){
            if (distr01(eng) > 0.5) child1[i] = (child1[i]-200)*(distr01(eng)*0.5+1) +200;
            else child1[i] = (child1[i]-200)*(1-distr01(eng)*0.5) +200;
        }
    }
    return child1;
}

vector<double> mutateAddTree(vector<double> child){
    int randomDepth = (int) (distr01(eng)*2)+1;
    int depthChild = find_depth(child);
    if (depthChild+randomDepth > 8) return child;
    while (randomDepth >= depthChild) {
        randomDepth = (int) (distr01(eng)*2)+1;
    }
    
    vector<double> randomTree = create_tree_full(randomDepth);
    int depth = depthChild-randomDepth+1;

    int r1 = (int) (distr01(eng)*(pow(2, depth-1))) + pow(2, depth-1);
    while (child[r1] == 99){
        r1 = (int) (distr01(eng)*(pow(2, depth-1))) + pow(2, depth-1);
    }
    int r2 = 1;
    int power = 0;
    int k = r2;
    for (int i = r1; i < 256; i=i*2){
        for (int l = 0; l<pow(2, power); l++){
            child[i+l] = randomTree[k+l];
        }
        k = k*2;
        power++;
    }
    return child;
}


vector<vector<double>> crossover(vector<double> parent1, vector<double> parent2){
    vector<double> child1 = parent1;
    vector<double> child2 = parent2;
    
    vector<vector<double>> res;
    if (parent1.size()==0 || parent2.size()==0){
        res.push_back(parent1);
        res.push_back(parent2);
        cout << "Parent size 0";
        return res;
    }
    int max = min( find_depth(parent1), find_depth(parent2));
    int r1 =(int) (distr01(eng)*(pow(2,max-1)-2))+2;
    int r2 = (int) (distr01(eng)*(pow(2,log2(r1+1))-2))+2;
    if (r2 == r1 && r2 != 1){
        r2 = (pow(2, (int)log2(r1)+1)-1) - (r1 - pow(2, (int)log2(r1)));
    }
    //cout << "r1: " << r1 << "; r2: " << r2 << endl;
    
    while (child1[r1] == 99 || child2[r2] == 99){
        r1 =(int) (distr01(eng)*(pow(2,max-1)-2))+2;
        r2 = (int) (distr01(eng)*(pow(2, (int) log2(r1)))) + pow(2, (int) log2(r1));
    }
    
    int power = 0;
    int k = r2;
    for (int i = r1; i < 256; i=i*2){
        for (int l = 0; l<pow(2, power); l++){
            swap(child1[i+l], child2[k+l]);
        }
        k = k*2;
        power++;
    }
    res.push_back(child1);
    res.push_back(child2);
    
    return res;
}

vector<vector<double>> crossover2(vector<double> parent1, vector<double> parent2){
    vector<double> child1 = parent1;
    vector<double> child2 = parent2;
    
    vector<vector<double>> res;
    if (parent1.size()==0 || parent2.size()==0){
        res.push_back(parent1);
        res.push_back(parent2);
        cout << "Parent size 0";
        return res;
    }
    int depth1 = find_depth(parent1);
    int depth2 = find_depth(parent2);
    int r1 =(int) (distr01(eng)*(pow(2,depth1-1)-2))+2;
    int r2 =(int) (distr02(eng)*(pow(2,depth2-1)-2))+2;
    
//    if (r2 == r1 && r2 != 1){
//        r2 = (pow(2, (int)log2(r1)+1)-1) - (r1 - pow(2, (int)log2(r1)));
//    }
        
    while (child1[r1] == 99 || child2[r2] == 99 || depth1-(int)log2(r1)+(int)log2(r2) > 8 ||depth2-(int)log2(r2)+(int)log2(r1) >8 ){
        //cout << "r1: " << r1 << "; r2: " << r2 << endl;
        r1 =(int) (distr01(eng)*(pow(2,depth1-1)-2))+2;
        r2 =(int) (distr02(eng)*(pow(2,depth2-1)-2))+2;
    }
    //cout << "r1: " << r1 << "; r2: " << r2 << endl;
    int power = 0;
    int k = r2;
    for (int i = r1; i < 256 && k < 256; i=i*2){
        for (int l = 0; l<pow(2, power); l++){
            child1[i+l] = parent2[k+l];
        }
        k = k*2;
        power++;
    }
    k = r1;
    power = 0;
    for (int i = r2; i < 256 && k < 256; i=i*2){
        for (int l = 0; l<pow(2, power); l++){
            child2[i+l] = parent1[k+l];
        }
        k = k*2;
        power++;
    }
    res.push_back(child1);
    res.push_back(child2);
    return res;
}

bool compareIndividuals(vector<double> ind1, vector<double> ind2){
    return calculate_fitness2(ind1) < calculate_fitness2(ind2);
}

void read_points(){
    ifstream input;
    input.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/points.txt");
    double err;
    while (input >> x >> y){
        Point a;
        a.x_coord = x;
        a.y_coord = y;
        //errorsHC.push_back(err);
        points.push_back(a);
    }
    input.close();
}

vector<double> create_tree_grow(int max_depth){
    vector<double> heap(256);
    size_t n = heap.size();
    
    for (int i = 0; i < n; i++){
        heap[i]  = 99;
    }
    for (int i = 1; i < pow(2,max_depth); i++){
        double current_node;
        if (int(log2(i)) < max_depth-1){
            current_node = 100+(int)(distr01(eng)*8);
            if (current_node == 100){
                current_node = 200 + (distr01(eng)*20-10);
            }
        }
        else {
            current_node = 100+(int) (distr01(eng)*2);
            if (current_node == 100){
                current_node = 200 + (distr01(eng)*20-10);
            }
        }
        if (i == 1){
            heap[i] = current_node;
        }
        else if (i%2==0){
            if (heap[(int)(i/2)] > 101){
                heap[i] = current_node;
            }
        }
        else if (i%2!=0){
            if (heap[(int)(i/2)] > 101){
                if (heap[(int)(i/2)] >= 106){
                    heap[i] = 99;
                }
                else if (heap[(int)(i/2)] <106){
                    heap[i] = current_node;
                }
            }
        }
    }
    return heap;
}

vector<double> create_tree_full(int max_depth){
    vector<double> heap(256);
    size_t n = heap.size();
    
    for (int i = 0; i < n; i++){
        heap[i]  = 99;
    }
    for (int i = 1; i < pow(2,max_depth); i++){
        double current_node;
        if (int(log2(i)) < max_depth-1){
            current_node = 100+(int)(distr01(eng)*6+2);
        }
        else {
            current_node = 100+(int) (distr01(eng)*2);
            if (current_node == 100){
                current_node = 200 + (distr01(eng)*20-10);
            }
        }
        if (i == 1){
            heap[i] = current_node;
        }
        else if (i%2==0){
            if (heap[(int)(i/2)] > 101){
                heap[i] = current_node;
            }
        }
        else if (i%2!=0){
            if (heap[(int)(i/2)] > 101){
                if (heap[(int)(i/2)] >= 106){
                    heap[i] = 99;
                }
                else if (heap[(int)(i/2)] <106){
                    heap[i] = current_node;
                }
            }
        }
    }
    return heap;
}

double calculate_function(double x, vector<double> heap){
    size_t n = heap.size();
    //if (n == 0) return 1000;
    vector<double> res(n);
    for (int i = 0; i < n; i++){
        res[i] = heap[i];
    }
    for (int i = n-1; i>=1; i--){
        switch ((int)heap[i])
        {
            case 99:
                break;
            case 101: res[i] = x;
                     break;
            case 102: res[i] = res[i*2]*res[i*2+1];
                    break;
            case 103: if (res[i*2+1] != 0) res[i] = res[i*2]/res[i*2+1];
                    else res[i] = 1;
                    break;
            case 104: res[i] = res[i*2]+res[i*2+1];
                    break;
            case 105: res[i] = res[i*2]-res[i*2+1];
                    break;
            case 106: res[i] = sin(res[i*2]);
                    break;
            case 107: res[i] = cos(res[i*2]);
                    break;
            default: res[i] = res[i]-200.0;
                    //cout << "const: " << res[i] << endl;
                     break;
        }
    }
    return res[1];
}

double f_x(double x){
    double f_x = 2*x*x + cos(x);
    return f_x;
}
void create_points(int n, int low, int high){
    ofstream output_points;
    output_points.open("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/points.txt");
    for (int i = 0; i < n-1; i++){
        double x = low + distr01(eng)*(high-low);
        double y = f_x(x);
        output_points << x << " " << y << endl;
    }
    output_points.close();
}
void plot_points(){
    Gnuplot gp("points");
    gp.set_xlabel("x");
    gp.set_ylabel("y");
    gp.set_title("Data points");
    gp.plotfile_xy("/Users/zhazira/Desktop/Evolution/Assignment2/Assignment2/points.txt", 1, 2, "Data points");
    gp.showonscreen();
    wait_for_key();
}
void plot_equation(){
    Gnuplot gp("lines");
    gp.set_xlabel("x");
    gp.set_ylabel("y");
    gp.set_title("Equation");
    gp.plot_equation("(sin(x))/x");
    gp.showonscreen();
    wait_for_key();
}
void wait_for_key ()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;
    
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}
