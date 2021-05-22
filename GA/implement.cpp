#include<stdio.h> 
#include<iostream> 
#include<stdlib.h> 
#include<unordered_map> 
#include<math.h> 
#include</Users/PremBhatia1/stdc.h> 
#include<random> 
#include<time.h> 
#define POOL 100   
#define ITERS 820             
#define BRK 90 

using namespace std ; 

double euclidean_dist(int x1, int y1, int x2, int y2) ; 

class Location
{
    public : 
        int X ; 
        int Y ; 
    
    public: 
        Location()
        {}

    public : 
        Location(int Xloc, int Yloc)
        {
            this->X = Xloc ; 
            this->Y = Yloc ; 
        }
 
        void assignLocation(int Xarg, int Yarg)
        {
            this->X = Xarg ; 
            this->Y = Yarg ; 
        } 

        double display(int mx, int my)
        {
            double D = euclidean_dist(X,Y,mx,my) ; 
            cout << "X_loc = " << X << " Y_loc = " << Y << " dist = " << D << endl  ; 
            return D ; 
        }
}; 

class Median
{
    public : 
        int X ; 
        int Y ; 
        int Pannels_to_median = 0 ; 
        double total_edge_weight = 0 ; 
        int desired_X_sum = 0 ; 
        int desired_Y_sum = 0 ; 
        unordered_map <int,Location> Panels ; /* List of Panels for that median mapped to an integer for O(1) access */
        std::vector<Location> panel_locations ; 
        unordered_map <int,double> weights_for_panels ; 

    public : 
        Median(int Xloc, int Yloc)
        {
            this->X = Xloc ; 
            this->Y = Yloc ;  
        }
 
        void assignPanel(int tag, Location loc, double Distance)
        {
            this->Panels[tag] = loc ; 
            this->weights_for_panels[tag] = Distance ; 
            this->panel_locations.push_back(loc) ; 
            this->Pannels_to_median++ ;  
            this->total_edge_weight = this->total_edge_weight + Distance ; 
            this->desired_X_sum = this->desired_X_sum + loc.X ; 
            this->desired_Y_sum = this->desired_Y_sum + loc.Y ; 
        }

        void display()
        {
            cout << "X_M = " << X <<  " Y_M = " << Y << " and Weight = " << this->total_edge_weight << " and panels = " << this->Pannels_to_median << endl ;   
        }

        int panelDisplay()
        {
            int x = 0 ; 
            cout << this->X << " :: " << this->Y << endl  ;
            double total = 0 ; 
            for(auto i:Panels)
            {
                x++ ; 
                //cout << i.first << "|" ; 
                //cout << this->X << " :: " << this->Y << endl  ; 
                total = total + i.second.display(X,Y) ; 
                cout << " ------- " << total << endl ; 
                //cout << endl ;   
            }
             
          return x ; 
        }
}; 

class Chromosome
{
    public : 
        std::vector<Median> mediansList ; 
        //int Panels_on_chrome[] ; 
        double totalEdgeWeight = 0 ; 

    public : 
        Chromosome()
        {}

    public : 
        Chromosome(std::vector<Median> List_of_Medians)
        {
            this->mediansList = List_of_Medians ; 
        }

    void display(int medians)
    {
        double cumulative = 0 ; 
        for(int kk = 0 ; kk < medians ; kk++)
        {
            cumulative = this->mediansList[kk].total_edge_weight + cumulative ; 
            this->mediansList[kk].display() ;  
            cout << " ---- " << cumulative << " ---- " << endl ;  
        }
        cout << endl ; 
    }
}; 

class IntPair
{ 
    public : 
        int first ; 
        int second ; 

    public : 
        IntPair()
        {}

    public : 
        IntPair(int f, int s)
    {
        this->first = f ; 
        this->second = s ; 
    }
}; 

int secondof(int a, int b, int c)
{
    if(a < b && a > c)
    {
        return a ; 
    }
    if(b > a && b < c)
    {
        return b ; 
    }
    return c ; 
}

IntPair* TS(int popSize, double prob)
{
    double P ; 
    int pairs[2] ;  
    srand(rand()) ;
    int Maxx = -1 ; 
    int Minn = 10000 ; 
    int Sec ; 
    for(int i = 0 ; i < 2 ; i++)
    {
        //Maxx = -1 ; 
        //Minn = 10000 ;
        //commenting this leads to generation of same chromes as parents which isnt much useful, yet it generates better results.  
        int first = rand()%popSize ;

        cout << "First = " << first << endl ; 

        if(first > Maxx)
        {
            Maxx = first ; 
        }

        if(first < Minn)
        {
            Minn = first ; 
        }

        int second = rand()%popSize ;

        cout << "Second = " << second << endl ; 

        if(second > Maxx)
        {
            Maxx = second ; 
        } 

        if(second < Minn)
        {
            Minn = second ; 
        }

        int third = rand()%popSize ;

        cout << "Third = " << third << endl ; 

        if(third > Maxx)
        {
            Maxx = third ; 
        }

        if(third < Minn)
        {
            Minn = third ; 
        }

        cout << "Min = " << Minn << endl ; 
        cout << "Maxx = " << Maxx << endl ; 
        cout << "Second = " << secondof(first,second,third) << endl ; 

        P = (double) rand()/RAND_MAX ; 

        cout << "P chosen is = " << P << endl ; 

        pairs[i] = Maxx ; 
        if(P >= prob)
        {
            pairs[i] = Minn ;  
        }
        if(P >= prob*(1-prob)  && P < prob)
        {
            cout << "P 2 Case" << endl ; 
            pairs[i] = secondof(first,second,third) ; 
        }
        
    }

    IntPair* pair = new IntPair(pairs[0], pairs[1]) ;  
    return pair ; 
}

double euclidean_dist(int x1, int y1, int x2, int y2)
{
    int ydiff = y2 - y1 ; 
    int xdiff = x2 - x1 ; 
    double dist = sqrt((double)(ydiff*ydiff + xdiff*xdiff)) ;  
    return dist ; 
}

void repairChrome(Chromosome* chromSoln, int P, int panels)
{
    cout << "Total OLD Weight = " << chromSoln->totalEdgeWeight << endl ; 
    double distPanels = 0 ; 
    std::vector<Median> ListMedians = chromSoln->mediansList ;
    cout << "P = " << P << endl ;  
    for(int k = 0 ; k < P ; k++)
    {
        double medianWeight = 0 ; 
        unordered_map<int,Location> panelsToMedian = chromSoln->mediansList.at(k).Panels ;
        cout << "k = " << k << endl ; 
        cout << "PREVIOUS - " << ListMedians.at(k).X << " - " << ListMedians.at(k).Y << endl ;  
        if(ListMedians.at(k).Pannels_to_median > 0)
        {
            ListMedians.at(k).X = ListMedians.at(k).desired_X_sum/ListMedians.at(k).Pannels_to_median ; 
            ListMedians.at(k).Y = ListMedians.at(k).desired_Y_sum/ListMedians.at(k).Pannels_to_median ; 
        } 
        cout << " - " << ListMedians.at(k).X << " - " << ListMedians.at(k).Y << endl ; 
        for(auto i : panelsToMedian)
        {
            distPanels = distPanels + euclidean_dist(i.second.X, i.second.Y, ListMedians.at(k).X, ListMedians.at(k).Y) ; 
            medianWeight = medianWeight + euclidean_dist(i.second.X, i.second.Y, ListMedians.at(k).X, ListMedians.at(k).Y) ; 
            //cout << "D = " << distPanels << endl ; 
        }
        ListMedians.at(k).total_edge_weight = medianWeight ;    
    }
    chromSoln->totalEdgeWeight = distPanels ; 
    chromSoln->mediansList = ListMedians ; 
    cout << "Total New Weight = " << chromSoln->totalEdgeWeight << endl ; 
}

Chromosome* createPopulation(Location* arrPanels, int medians, int panels, int X_Size, int Y_Size , int PanelLocs[]) 
{
    Chromosome* chromes = new Chromosome[POOL] ; 
    int medianEntries[X_Size*Y_Size] ; 
    srand(time(0)) ; 
    for(int j = 0 ; j < POOL ; j++)
    {  
        std::vector<Median> medians_for_pool ;
        for(int l = 0 ; l < medians ; l++)
        {
            int Xm = rand()%X_Size ; int Ym = rand()%Y_Size ; 
            Median median_point(Xm,Ym) ;
            medians_for_pool.push_back(median_point) ; 
            medianEntries[Xm*Y_Size + Ym] = 1 ;  
        }
        for(int p = 0 ; p < panels ; p++)
        {
            int MED = -1 ; 
            double DIST = -1 ; 
            double min_dist = std::numeric_limits<double>::max() ;     
            for(int m = 0 ; m < medians ; m++) 
            { 
                Median med = medians_for_pool.at(m) ;  
                double distance = euclidean_dist(arrPanels[p].X, arrPanels[p].Y, med.X, med.Y) ;
                if(distance < min_dist)
                {
                    min_dist = distance ;  
                    MED = m ;  
                    DIST = distance ; 
                } 
            }  
            medians_for_pool.at(MED).assignPanel(arrPanels[p].X*Y_Size + arrPanels[p].Y, arrPanels[p],DIST) ;  
            chromes[j].totalEdgeWeight = chromes[j].totalEdgeWeight + DIST ; 
            cout << "Panel = P " << p << " assigned to Median = " << MED << " Chrome["<<j<<"] "<< DIST << endl ;  
        }
        chromes[j].mediansList = medians_for_pool ; 
        repairChrome(chromes+j, medians, panels) ; 
    } 
    return chromes ; 
}

void swap(Chromosome* Population, int first, int second)
{
    Chromosome Temp ; 
    Temp = Population[first] ; 
    Population[first] = Population[second] ; 
    Population[second] = Temp ; 
}

IntPair* select_parent(Chromosome* GA_population)
{
    double cdf = 0 ; 
    int pairs[2] ; 
    pairs[0] = 0 ; 
    pairs[1] = 1 ; 
    double parents[2] ; 
    parents[0] = (double) rand()/RAND_MAX ; 
    //parents[0] = parents[0]/2 ; 
    //parents[1] = (double) rand()/RAND_MAX ; 
    //parents[1] = parents[0] + 0.5 ; 
    parents[1] = (double) rand()/RAND_MAX ; 
    double total = 0 ; 

    for(int i = 0 ; i < POOL ; i++)
    {
        total = total + GA_population[i].totalEdgeWeight ; 
    }

    cout << "For Selection !" << endl ; 

    for(int j = POOL ; j > 0 ; j--)
    {
        cdf = cdf + GA_population[j-1].totalEdgeWeight ; 
        //if(parents[0] <  cdf/total*0.75 )
        if(parents[0] <  cdf/total) 
        {
            pairs[0] = j - 1 ; 
            cout << " - " << parents[0] << endl ;
            break ; 
	    cout << " Decider = " << cdf/total*0.75 << " where cdf = " << cdf << " and total = " << total << endl ;   
            cout << "Assigned pairs[0] " << j - 1 << " as -->>>>> " << GA_population[j-1].totalEdgeWeight << endl ;  
        }
    }
    cdf = 0 ; 
    for(int j = POOL ; j > 0 ; j--)
    {
        cdf = cdf + GA_population[j-1].totalEdgeWeight ; 
        //if(parents[1] < cdf/total*0.75 )
        if(parents[1] < cdf/total) 
        {
            pairs[1] = j - 1 ;
            if(pairs[1] == 0)
            {
                pairs[1] = 1 ; 
            }  
            cout << " - " << parents[1] << endl ;  
            break ; 
	    cout << " Decider = " << cdf/total*0.75 << " where cdf = " << cdf << " and total = " << total << endl ;  
            cout << "Assigned pairs[1] " << j - 1 << " as -->>>>> " << GA_population[j-1].totalEdgeWeight << endl ;  
        }
    }

    IntPair* generators = new IntPair(pairs[0],pairs[1]) ;  
    return generators ; 
}

Chromosome* allocatePanels(std::vector<Median> median_list_for_allocation, Location* arrPanels, int panels, int medians, int X_Size, int Y_Size)
{
    Chromosome* OffSpring = new Chromosome() ; 
    for(int i = 0 ; i < panels ; i++)
    {
        int MED = -1 ; 
        double DIST = -1 ; 
        double min_dist = std::numeric_limits<double>::max() ; 
        for(int m = 0 ; m < medians ; m++)
        {
            Median mFromChild = median_list_for_allocation.at(m) ; 
            double distance = euclidean_dist(arrPanels[i].X, arrPanels[i].Y, mFromChild.X, mFromChild.Y) ; 
            if(distance < min_dist)
            {
                min_dist = distance ; 
                MED = m ; 
                DIST = distance ; 
            }
        }
        median_list_for_allocation.at(MED).assignPanel(arrPanels[i].X*Y_Size + arrPanels[i].Y, arrPanels[i], DIST) ;
        OffSpring->totalEdgeWeight = OffSpring->totalEdgeWeight + DIST ;  
    }
    OffSpring->mediansList = median_list_for_allocation ; 
    return OffSpring ; 
}

Chromosome* crossover(std::vector<Median> Mlist_first, std::vector<Median> Mlist_second, int P, Location* arrPanels, int panels, int X_Size, int Y_Size)
{
    std::vector<Median> offspring_median_list ; 
    srand(rand()) ;  
    int hinge = rand()%P ; 
    cout << "hinged at -> " << hinge << endl ; 
    for(int start = 0 ; start < hinge ; start++)
    {
        Median med_for_cross(Mlist_first.at(start).X, Mlist_first.at(start).Y) ; 
        offspring_median_list.push_back(med_for_cross) ; 
    }
    for(int start = hinge ; start < P ; start++)
    {
        Median med_for_cross(Mlist_second.at(start).X, Mlist_second.at(start).Y) ; 
        offspring_median_list.push_back(med_for_cross) ; 
    }
    cout << "HINGE = " << hinge << endl ; 
    Chromosome* offSpring = allocatePanels(offspring_median_list, arrPanels, panels, P, X_Size, Y_Size) ; 
    return offSpring ; 
}

int partition (Chromosome* chromes, int low, int high)
{
    int pivot = chromes[high].totalEdgeWeight ; // pivot
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
 
    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (chromes[j].totalEdgeWeight < pivot)
        {
            i++; // increment index of smaller element
            swap(chromes, i, j); 
        }
    }
    swap(chromes, i+1, high); 
    return (i + 1);
}
 
/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort(Chromosome* chromes, int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
        at right place */
        int pi = partition(chromes, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(chromes, low, pi - 1);
        quickSort(chromes, pi + 1, high);
    }
}

void Sort_Chromes(Chromosome* chromes)
{
    //quickSort(chromes, 0, POOL-1) ;   
    for(int j = 0 ; j < POOL ; j++)
    {
        for(int i = 0 ; i < POOL - j - 1 ; i++)
        {
            if(chromes[i].totalEdgeWeight > chromes[i+1].totalEdgeWeight) 
            {
                swap(chromes, i, i+1) ; 
            }
        }
    } 
}

void Check_Order(Chromosome* chromes)
{
    int check = 1 ; 
    for(int i = 0 ; i < POOL-1 ; i++)
    {
        if(chromes[i].totalEdgeWeight > chromes[i+1].totalEdgeWeight)
        check = 0 ; 
    }
    if(check == 1)
    {
        cout << " - - - SORTED - - -" << endl ;  
    }
}


IntPair* TS_2(int Pool,int Tsize)
{
    int First = 0 ; 
    int Min = 100000 ; 
    for(int i = 0 ; i < Tsize ; i++)
    {
        int random = rand()%Pool ;
        cout << "       -> " << random << endl ; 
        if(Min > random)
        {
            Min = random ;
            First = Min ;  
        }
    }
    int Second = 1 ; 
    Min = 100000 ;
    for(int i = 0 ; i < Tsize ; i++)
    {
        int random = rand()%Pool ;
        cout << "       -> " << random << endl ; 
        if(Min > random)
        {
            Min = random ;
            Second = Min ;  
        }
    }
    cout << "First = " << First << " and " << " and Second = " << Second << endl ;  
    IntPair* gen = new IntPair(First,Second) ; 
    return gen ; 
}

void Start_the_GA(double Best, Chromosome* GA_population, int P, Location* Panels, int N, int X, int Y)
{
    ofstream myfile("example.txt",ios::app) ; 
    int converge = 0 ; 
    srand(time(0)) ; 
    for(int iters = 0 ; iters < ITERS ; iters++)
    {
        cout << "iteration number is -> " << iters << endl ; 
        if(iters > ITERS-5) 
        {
            myfile << "Iterations : " << iters << " Solution : " << GA_population[0].totalEdgeWeight << endl ; 
            break ; 
	    }
        Sort_Chromes(GA_population) ; 
        // This function is just for checking if sorted or not 
        //Check_Order(GA_population); 
        if(GA_population[0].totalEdgeWeight != Best)
        {
            cout << "Not Equal as : best = " << Best << " and GA[0] = " << GA_population[0].totalEdgeWeight << endl ; 
            cout << Best - GA_population[0].totalEdgeWeight << endl ; 
            converge = 0 ; 
            Best = GA_population[0].totalEdgeWeight ; 
        }
        else{
            cout << Best - GA_population[0].totalEdgeWeight << endl ; 
            cout << "Check Here" << endl ; 
	    converge++ ; 
        }

        if(converge > 10)
        {
            cout << " ------------ CONVERGENCE " <<  converge << " STARTED AS MAX --> " << GA_population[0].totalEdgeWeight << endl ; 
            if(converge > BRK) 
            {
                myfile << "Iterations : " << iters << " Solution : " << GA_population[0].totalEdgeWeight << endl ; 
                //myfile << "case     " << N << "     " << P << "     " << GA_population[0].totalEdgeWeight << endl ;  
                break ; 
            }
        }


        //IntPair* cross = TS(POOL, 0.05) ; 
        /*This is 3 Way Tournament */
        IntPair* cross = TS_2(POOL,3) ;  

        /*This code makes use of SUS or FPS */
        //IntPair* parents = select_parent(GA_population) ;
        /* SUS */
        //IntPair* cross = parents ;  

        cout << "C -> first = " << cross->first << endl ; 
        cout << "C -> second = " << cross->second << endl ; 

        cout << "P_1 weight = " << GA_population[cross->first].totalEdgeWeight << endl ; 
        cout << "P_2 weight = " << GA_population[cross->second].totalEdgeWeight << endl ; 

        Chromosome* off_spring = crossover(GA_population[cross->first].mediansList, GA_population[cross->second].mediansList, P, Panels, N, X, Y) ;
        cout << "Crossover new offspring Weight = " << off_spring->totalEdgeWeight << endl ; 
        GA_population[POOL-1] = *off_spring ;  
        int p_mutation = rand()%100 ;  
        cout << "-  -   -   -   -   " << p_mutation << endl ; 
        if(p_mutation > 20)
        {
            cout << "+  +   +   + So Enters" << endl ; 
            int index = rand()%P ; 
            repairChrome(GA_population+index, P, N) ; 
        }
    }
}

int main(int argc, char** argv)
{
    clock_t start, end ; 
    start = clock() ; 
    cout << "Test" << endl ;
    int X ; int Y ; int N ; int P ; int panelsX ; int panelsY ;  int counts ; /* X_length, Y_length, No. of Panels, Medians */
    cin >> X ; 
    cin >> Y ;
    cin >> N ;  
    cin >> P ;   
    // Array for all panels . Datatype is location. 
    Location* Panels = new Location[N] ;
    // 2D array for all Panel locations. 
    cout << "Array mem. allocation" << endl ; 
    int* PanelLocations = new int[X*Y] ;   
    cout << "Not reached for large values of X,Y" << endl ; 
    int obstacles[N] ; 
    cout << "1st Loop starting" << endl ; 
    for(int ii = 0 ; ii < N ; ii++)
    {
        cin >> panelsX ; 
        cin >> panelsY ;
        if(rand()%100 < 2)
        {
            cout << "       -- Checked It From Here     !! " << ii << endl ;  
            panelsX = int(panelsX*0.036 + panelsX)%100 ; 
            panelsY = int(panelsY*0.036 + panelsY)%100 ; 
            cout << panelsX << " -- " << panelsY << endl ; 
        } 
        PanelLocations[panelsX*Y + panelsY] = 1 ; 
        Panels[ii] = Location(panelsX,panelsY) ; 
    }  
    Chromosome* population = createPopulation(Panels,P,N,X,Y,PanelLocations) ;
    cout << "Created Population " << endl ; 

    Sort_Chromes(population) ; 

    cout << "Done Sorting " << endl ;

    for(int j = 0 ; j < POOL ; j++)
    {
        cout << j << ". " << population[j].totalEdgeWeight << endl ; 
    }
    cout << endl ; 

    population[0].display(P) ; 
    population[0].mediansList.at(0).panelDisplay() ; 
    population[0].mediansList.at(P-1).panelDisplay() ;   

    double Best = 0 ; 

    Start_the_GA(Best, population, P, Panels, N, X, Y) ;  

    Sort_Chromes(population) ;

    for(int j = 0 ; j < POOL ; j++)
    {
        cout << j << ". " << population[j].totalEdgeWeight << endl ;
    } 
 

    population[0].display(P) ;
    int length_m = population[0].mediansList.size() ; 
    population[0].mediansList.at(0).panelDisplay() ;  
    population[0].mediansList.at(1).panelDisplay() ;
    population[0].mediansList.at(2).panelDisplay() ;  
    population[0].mediansList.at(3).panelDisplay() ;
    population[0].mediansList.at(4).panelDisplay() ;  

    cout << "Total number of med indices = " << length_m << endl ;   

    /* 
    int total_meds = 0 ; 
    double for_check = 0 ; 
    for(int m = 0 ; m < P ; m++)
    {
        total_meds = total_meds + population[8].mediansList.at(m).panelDisplay() ;   
        cout << total_meds ;
        for_check = for_check + population[8].mediansList.at(m).total_edge_weight ; 
        cout << " " << population[8].mediansList.at(m).total_edge_weight << " - " << for_check << endl ;  
    } 
    cout << "T : " << population[8].totalEdgeWeight << endl ; 
    cout << total_meds << endl ; 

    //swap(population, 2, 8) ; 
    cout << "AFTER SWAP -> " << endl ; 
    total_meds = 0 ; 
    for_check = 0 ; 
    for(int m = 0 ; m < P ; m++)
    {
        total_meds = total_meds + population[2].mediansList.at(m).panelDisplay() ;   
        cout << total_meds ;
        for_check = for_check + population[2].mediansList.at(m).total_edge_weight ; 
        cout << " " << population[2].mediansList.at(m).total_edge_weight << " - " << for_check << endl ;  
    } 
    cout << "T : " << population[2].totalEdgeWeight << endl ; 
    cout << total_meds << endl ; 
    cout << endl ;  

    Chromosome* offSpring_2_8 = crossover(population[8].mediansList, population[2].mediansList, P, Panels, N, X, Y) ;

    total_meds = 0 ; 
    for_check = 0 ; 
    for(int m = 0 ; m < P ; m++)
    {
        total_meds = total_meds + offSpring_2_8->mediansList.at(m).panelDisplay() ;   
        cout << total_meds ; 
        for_check = for_check + offSpring_2_8->mediansList.at(m).total_edge_weight ; 
        cout << " " << offSpring_2_8 ->mediansList.at(m).total_edge_weight << " - " << for_check << endl ;  
    }
    cout << "T : " << offSpring_2_8->totalEdgeWeight << endl ; 
    cout << total_meds << endl ; 
    cout << endl ; 

    for(int j = 0 ; j < POOL ; j++)
    {
        cout << j << ". " << population[j].totalEdgeWeight << endl ;
        population[j].display(P) ;  
    }
    cout << endl ; 

    Sort_Chromes(population) ; 

    for(int j = 0 ; j < POOL ; j++)
    {
        cout << j << ". " << population[j].totalEdgeWeight << endl ;
        population[j].display(P) ;  
    }
    */ 
    delete[] Panels ;  
    delete[] population ;  
    cout << "Yes -> The GA - : Working Fine !! :) " << endl ; 
    end = clock() ; 
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
    ofstream Timefile("time_info.txt",ios::app) ;
    Timefile << "Time taken by program for N = " << N << " and P = " << P << " is : " << fixed << time_taken << setprecision(2) ; 
    Timefile << " sec " << endl ;     
    cout << "Done - GASPV_soln_mpj" << endl ; 
    vector<int> vec ; 
    return 0 ; 
}
