// ConsoleApplication2.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
/*
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <random>

class Agent {
public:
    char status;
    int timeElapsed;
    int dE, dI, dR;

    Agent(char s, int de, int di, int dr) : status(s), dE(de), dI(di), dR(dr), timeElapsed(0) {}

    void updateStatus() {
        if (status == 'E' && timeElapsed >= dE) {
            status = 'I';
            timeElapsed = 0;
        }
        else if (status == 'I' && timeElapsed >= dI) {
            status = 'R';
            timeElapsed = 0;
        }
        else if (status == 'R' && timeElapsed >= dR) {
            status = 'S';
            timeElapsed = 0;
        }
        else {
            ++timeElapsed;
        }
    }
};

class Simulation {
public:
    std::vector<std::vector<Agent>> grid;
    int totalDays;
    int gridSize;
    int maxInfections;

    Simulation(int size, int days, int maxInf)
        : gridSize(size), totalDays(days), maxInfections(maxInf) {
        grid.resize(gridSize, std::vector<Agent>(gridSize, Agent('S', 3, 5, 7)));

        // Initialisation aléatoire des agents
        initializeAgents();
    }

    void initializeAgents() {
        // Générateur de nombres aléatoires (Mersenne Twister)
        std::mt19937 rng(static_cast<unsigned>(time(nullptr)));

        // Distribution exponentielle avec les paramètres 3, 7, et 365
        std::exponential_distribution<double> exp_dE(1.0 / 3.0);
        std::exponential_distribution<double> exp_dI(1.0 / 7.0);
        std::exponential_distribution<double> exp_dR(1.0 / 365.0);

        int infectedCount = 0;
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                double randomValue = rng() / static_cast<double>(rng.max());

                if (infectedCount < maxInfections) {
                    if (randomValue < 0.0001) { // 0.01% chance of being infected
                        grid[i][j] = Agent('I', exp_dE(rng), exp_dI(rng), exp_dR(rng));
                        ++infectedCount;
                    }
                }
                else {
                    if (randomValue < 0.999) { // 99.9% chance of being susceptible
                        grid[i][j] = Agent('S', exp_dE(rng), exp_dI(rng), exp_dR(rng));
                    }
                }
            }
        }
    }
     void run() {
        for(int k=0;k<2;k++){
            std::stringstream filename;
            filename << "results_" << k << ".CSV";
            std::ofstream file(filename.str());
            if (file.is_open()) {
                file << "Day          S          E          I          R" << std::endl;
                for (int day = 0; day < totalDays; ++day) {
                    for (int i = 0; i < gridSize; ++i) {
                        for (int j = 0; j < gridSize; ++j) {
                            Agent& agent = grid[i][j];
                            agent.updateStatus();
                            // Propagation de la maladie - Ajoutez votre logique ici
                        }
                    }
                    // Propagation de la maladie
                    for (int i = 0; i < gridSize; ++i) {
                        for (int j = 0; j < gridSize; ++j) {
                            Agent& agent = grid[i][j];
                            if (agent.status == 'S') {
                                int infectedCount = 0;

                                // Comptez le nombre d'agents infectieux dans le voisinage (8 cases autour)
                                for (int dx = -1; dx <= 1; ++dx) {
                                    for (int dy = -1; dy <= 1; ++dy) {
                                        int x = (i + dx + gridSize) % gridSize;
                                        int y = (j + dy + gridSize) % gridSize;
                                        if (grid[x][y].status == 'I') {
                                            ++infectedCount;
                                        }
                                    }
                                }

                                // Calculez la probabilité d'infection
                                double p = 1.0 - exp(-0.5 * infectedCount);
                                if (rand() / double(RAND_MAX) < p) {
                                    agent.status = 'E';
                                    agent.timeElapsed = 0;
                                }
                            }
                        }
                    }
                    int sCount = 0, eCount = 0, iCount = 0, rCount = 0;
                    for (int i = 0; i < gridSize; ++i) {
                        for (int j = 0; j < gridSize; ++j) {
                            Agent& agent = grid[i][j];
                            if (agent.status == 'S') ++sCount;
                            else if (agent.status == 'E') ++eCount;
                            else if (agent.status == 'I') ++iCount;
                            else if (agent.status == 'R') ++rCount;
                        }
                    }
                    file << day << "          " << sCount << "          " << eCount << "          " << iCount << "          " << rCount << std::endl;
                }
                file.close();
            }
        }
    }

  
};

int main() {
    Simulation simulation(300, 730, 20);
    simulation.run();
   // simulation.saveResults();
    return 0;
}*/
/*
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <random>

class Agent {
public:
    char status;
    int timeElapsed;
    int dE, dI, dR;

    Agent(char s, int de, int di, int dr) : status(s), dE(de), dI(di), dR(dr), timeElapsed(0) {}

    void updateStatus() {
        if (status == 'E' && timeElapsed >= dE) {
            status = 'I';
            timeElapsed = 0;
        }
        else if (status == 'I' && timeElapsed >= dI) {
            status = 'R';
            timeElapsed = 0;
        }
        else if (status == 'R' && timeElapsed >= dR) {
            status = 'S';
            timeElapsed = 0;
        }
        else {
            ++timeElapsed;
        }
    }
};

class Simulation {
public:
    std::vector<std::vector<Agent>> grid;
    int totalDays;
    int gridSize;
    int maxInfections;
    int totalAgents;

    Simulation(int size, int days, int maxInf, int total) : gridSize(size), totalDays(days), maxInfections(maxInf), totalAgents(total) {
        grid.resize(gridSize, std::vector<Agent>(gridSize, Agent('S', 3, 5, 7)));

        // Initialisation des agents
        initializeAgents();
    }*/

    /*void initializeAgents() {
        // Générateur de nombres aléatoires (Mersenne Twister)
        std::mt19937 rng(static_cast<unsigned>(time(nullptr)));

        // Distribution exponentielle avec les paramètres 3, 7, et 365
        std::exponential_distribution<double> exp_dE(1.0 / 3.0);
        std::exponential_distribution<double> exp_dI(1.0 / 7.0);
        std::exponential_distribution<double> exp_dR(1.0 / 365.0);

        int infectedCount = 0;
        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                double randomValue = rng() / static_cast<double>(rng.max());

                if (infectedCount < maxInfections) {
                    if (randomValue < static_cast<double>(maxInfections) / totalAgents) {
                        grid[i][j] = Agent('I', exp_dE(rng), exp_dI(rng), exp_dR(rng));
                        ++infectedCount;
                    }
                }
                else {
                    if (randomValue < static_cast<double>(totalAgents - maxInfections) / totalAgents) {
                        grid[i][j] = Agent('S', exp_dE(rng), exp_dI(rng), exp_dR(rng));
                    }
                }
            }
        }
    }
    void initializeAgents() {
        // Générateur de nombres aléatoires (Mersenne Twister)
        std::mt19937 rng(static_cast<unsigned>(time(nullptr)));

        // Distribution exponentielle avec les paramètres 3, 7, et 365
        std::exponential_distribution<double> exp_dE(1.0 / 3.0);
        std::exponential_distribution<double> exp_dI(1.0 / 7.0);
        std::exponential_distribution<double> exp_dR(1.0 / 365.0);

        std::vector<std::pair<int, int>> infectedLocations;

        // Générez 20 emplacements uniques pour les individus infectés
        while (infectedLocations.size() < maxInfections) {
            int x = rng() % gridSize;
            int y = rng() % gridSize;
            if (std::find(infectedLocations.begin(), infectedLocations.end(), std::make_pair(x, y)) == infectedLocations.end()) {
                infectedLocations.push_back(std::make_pair(x, y));
            }
        }

        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                if (std::find(infectedLocations.begin(), infectedLocations.end(), std::make_pair(i, j)) != infectedLocations.end()) {
                    grid[i][j] = Agent('I', exp_dE(rng), exp_dI(rng), exp_dR(rng));
                }
                else {
                    grid[i][j] = Agent('S', exp_dE(rng), exp_dI(rng), exp_dR(rng));
                }
            }
        }
    }*/
    /*void initializeAgents() {
        // Générateur de nombres aléatoires (Mersenne Twister)
        std::mt19937 rng(static_cast<unsigned>(time(nullptr)));

        // Distribution exponentielle avec les paramètres 3, 7, et 365
        std::exponential_distribution<double> exp_dE(1.0 / 3.0);
        std::exponential_distribution<double> exp_dI(1.0 / 7.0);
        std::exponential_distribution<double> exp_dR(1.0 / 365.0);

        // Liste des emplacements pour les individus infectés
        std::vector<std::pair<int, int>> infectedLocations;

        // Répartissez aléatoirement les 20 individus infectés
        while (infectedLocations.size() < 20) {
            int x = rng() % gridSize;
            int y = rng() % gridSize;
            if (std::find(infectedLocations.begin(), infectedLocations.end(), std::make_pair(x, y)) == infectedLocations.end()) {
                infectedLocations.push_back(std::make_pair(x, y));
                grid[x][y] = Agent('I', exp_dE(rng), exp_dI(rng), exp_dR(rng));
            }
        }

        // Répartissez aléatoirement les 19980 individus susceptibles
        for (int i = 0; i < 19980; ++i) {
            int x, y;
            do {
                x = rng() % gridSize;
                y = rng() % gridSize;
            } while (grid[x][y].status != 'S');
            grid[x][y] = Agent('S', exp_dE(rng), exp_dI(rng), exp_dR(rng));
        }
    }



    void run() {
        for (int k = 0; k < 2; k++) {
            std::stringstream filename;
            filename << "results_" << k << ".csv";
            std::ofstream file(filename.str());
            if (file.is_open()) {
                file << "Day, S, E, I, R" << std::endl;
                for (int day = 0; day < totalDays; ++day) {
                    for (int i = 0; i < gridSize; ++i) {
                        for (int j = 0; j < gridSize; ++j) {
                            Agent& agent = grid[i][j];
                            agent.updateStatus();
                        }
                    }
                    for (int i = 0; i < gridSize; ++i) {
                        for (int j = 0; j < gridSize; ++j) {
                            Agent& agent = grid[i][j];
                            if (agent.status == 'S') {
                                int infectedCount = 0;
                                for (int dx = -1; dx <= 1; ++dx) {
                                    for (int dy = -1; dy <= 1; ++dy) {
                                        int x = (i + dx + gridSize) % gridSize;
                                        int y = (j + dy + gridSize) % gridSize;
                                        if (grid[x][y].status == 'I') {
                                            ++infectedCount;
                                        }
                                    }
                                }
                                double p = 1.0 - exp(-0.5 * infectedCount);
                                if (rand() / static_cast<double>(RAND_MAX) < p) {
                                    agent.status = 'E';
                                    agent.timeElapsed = 0;
                                }
                            }
                        }
                    }
                    int sCount = 0, eCount = 0, iCount = 0, rCount = 0;
                    for (int i = 0; i < gridSize; ++i) {
                        for (int j = 0; j < gridSize; ++j) {
                            Agent& agent = grid[i][j];
                            if (agent.status == 'S') ++sCount;
                            else if (agent.status == 'E') ++eCount;
                            else if (agent.status == 'I') ++iCount;
                            else if (agent.status == 'R') ++rCount;
                        }
                    }
                    file << day << "  " << sCount << "  " << eCount << "  " << iCount << "  " << rCount << std::endl;
                }
                file.close();
            }
        }
    }
};

int main() {
    int gridSize = 300;
    int totalAgents = 20000; // Nombre total d'individus
    int totalDays = 730;
    int maxInfections = 20;

    // Assurez-vous que totalAgents ne dépasse pas la taille de la grille
    if (totalAgents > gridSize * gridSize) {
        std::cerr << "Le nombre total d'individus dépasse la taille de la grille. Réduisez le nombre d'individus ou agrandissez la grille." << std::endl;
        return 1;
    }

    Simulation simulation(gridSize, totalDays, maxInfections, totalAgents);
    simulation.run();

    return 0;
}*/
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <string>

const int TAILLE_GRILLE = 300;
const int NB_INDIVIDUS = 20000;
const int NB_INFECTES_INIT = 20;

double negExp(double inMean, std::mt19937& gen) {
    std::exponential_distribution<> distrib(1.0 / inMean);
    return distrib(gen);
}

enum Statut { S, E, I, R };

class Individu {
public:
    Statut statut;
    int tempsStatut;
    int dE, dI, dR;
    int posX, posY;

    Individu() : statut(S), tempsStatut(0), dE(0), dI(0), dR(0), posX(0), posY(0) {}

    void deplacer() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, TAILLE_GRILLE - 1);

        posX = distrib(gen);
        posY = distrib(gen);
    }

    void miseAJourStatut() {
        switch (statut) {
        case E:
            if (++tempsStatut > dE) { statut = I; tempsStatut = 0; }
            break;
        case I:
            if (++tempsStatut > dI) { statut = R; tempsStatut = 0; }
            break;
        case R:
            if (++tempsStatut > dR) { statut = S; tempsStatut = 0; }
            break;
        default:
            break;
        }
    }
};

class Grille {
private:
    std::vector<Individu> individus;
    std::vector<std::vector<std::vector<int>>> cellules;
    std::ofstream file;

public:
    Grille() : cellules(TAILLE_GRILLE, std::vector<std::vector<int>>(TAILLE_GRILLE)) {}

    void initialiser(int num_simulation) {
        file.open("simulation_" + std::to_string(num_simulation) + ".csv");
        file << "S,E,I,R\n";

        std::random_device rd;
        std::mt19937 gen(rd());

        for (int i = 0; i < NB_INDIVIDUS; ++i) {
            Individu ind;
            ind.deplacer();

            ind.dE = static_cast<int>(negExp(3.0, gen));
            ind.dI = static_cast<int>(negExp(7.0, gen));
            ind.dR = static_cast<int>(negExp(365.0, gen));

            if (i < NB_INFECTES_INIT) {
                ind.statut = I;
            }

            individus.push_back(ind);
            cellules[ind.posX][ind.posY].push_back(i);
        }
    }

    void iteration(int jour) {
        for (int i = 0; i < individus.size(); ++i) {
            auto& ind = individus[i];
            auto& ancienneCellule = cellules[ind.posX][ind.posY];
            ancienneCellule.erase(std::remove(ancienneCellule.begin(), ancienneCellule.end(), i), ancienneCellule.end());

            ind.deplacer();
            cellules[ind.posX][ind.posY].push_back(i);

            ind.miseAJourStatut();
        }

        for (auto& ind : individus) {
            if (ind.statut == S) {
                int ni = compterInfectesVoisinage(ind.posX, ind.posY);
                double proba = 1 - exp(-0.5 * ni);

                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> distrib(0.0, 1.0);
                if (distrib(gen) < proba) {
                    ind.statut = E;
                }
            }
        }

        int nbS = 0, nbE = 0, nbI = 0, nbR = 0;
        for (const auto& ind : individus) {
            switch (ind.statut) {
            case S: nbS++; break;
            case E: nbE++; break;
            case I: nbI++; break;
            case R: nbR++; break;
            }
        }
        file  << nbS << "," << nbE << "," << nbI << "," << nbR << "\n";
    }

    int compterInfectesVoisinage(int x, int y) {
        int count = 0;
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int nx = (x + dx + TAILLE_GRILLE) % TAILLE_GRILLE;
                int ny = (y + dy + TAILLE_GRILLE) % TAILLE_GRILLE;

                for (auto idx : cellules[nx][ny]) {
                    if (individus[idx].statut == I) {
                        count++;
                    }
                }
            }
        }
        return count;
    }

    void fermerFichier() {
        file.close();
    }
};

int main() {
    for (int sim = 1; sim <= 100; sim++) {
        Grille grille;
        grille.initialiser(sim);

        for (int jour = 0; jour < 730; ++jour) {
            grille.iteration(jour);
        }

        grille.fermerFichier();
    }

    return 0;
}
