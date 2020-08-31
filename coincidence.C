#include "transformations.h"

//returns largest element in int array
double maxElement(const double* arr, const int& size)
{
    double max = arr[0];

    for (int i = 0; i < size; ++i)
    {
        if(arr[i] > max) max = arr[i];
    }

    return max;
}

//templated struct for two arrays of same size used in sortSplitted
template<int N>
struct XYarr
{
    double x[N],y[N];
    int size() {return N;}

    //manually copying arrays in constructor
    XYarr(const double (&x_)[N], const double (&y_)[N]) 
    {
        memcpy(&x,x_,sizeof (double)*N);
        memcpy(&y,y_,sizeof (double)*N);
    }
};

//base event class
struct Event
{
    char det;           //detector code
    int id;             //ID of the event
    double mjd;         //modified julian date
    double ra, dec;     //right ascension and declination in degrees
    double resolution;  //circular resolution in degrees
};

//structure for loading data from IceCube catalog of alerts
struct IceCubeEvent : Event
{
    //unknown errors will be set to zero
    double ra_err_high, ra_err_low;     //right ascension 90% confidence interval (asymetric error size)
    double dec_err_high, dec_err_low;  //declination asymetric error size

    IceCubeEvent(char det, int id, double mjd, double ra, double dec, double reh, double rel,
                 double deh, double del) : Event{det, id, mjd, ra, dec, std::nan("")}, ra_err_high(reh),
                 ra_err_low(rel), dec_err_high(deh), dec_err_low(del) {}
};

//structure for loading general data from IceCube
struct IceCubeEvent2 : Event
{
    double logEnergy;

    IceCubeEvent2(char det, int id, double mjd, double ra, double dec, double res, double lE) 
        : Event{det, id, mjd, ra, dec, res}, logEnergy(lE) {}
};

//structure for loading data from ANTARES catalog
struct AntaresEvent : Event
{
    int nhits; //detected number of hits

    AntaresEvent(char det, int id, double mjd, double ra, double dec, int nhits, double res) 
        : Event{det, id, mjd, ra, dec, res}, nhits(nhits) {}
};

//structure for loading data from Baikal-GVD
struct BaikalEvent : Event
{
    int nhits;
    double energy;
    double energySigma;

    BaikalEvent(char det, int id, double mjd, double ra, double dec, double res, int nhits, 
                double e, double eS) : Event{det, id, mjd, ra, dec, res}, nhits(nhits), 
                energy(e), energySigma(eS) {}
};

//output structure with coincidence parameters
struct Coincidence
{
    std::string det;            //detector codes concatenated
    int id, id2;                //IDs of coinciding events
    double timeDev, angleDev;   //time and angle difference
};

//operator overloading for appending std::vector
template<typename T>
void operator+=(std::vector<T>& v1, const std::vector<T>& v2)
{
    v1.insert(v1.end(), v2.begin(), v2.end());
}

//operator overloading for printing
std::ostream& operator<<(std::ostream& stream, const Event& ev)
{
    stream << "Detector code: " << ev.det << "\nID: " << ev.id << '\n';
    stream << "MJD: " << ev.mjd << " RA: " << ev.ra << " Dec: " << ev.dec << '\n';
    stream << "Resolution: " << ev.resolution << '\n';

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const IceCubeEvent& ev)
{
    stream << "ID: " << ev.id << '\n';
    stream << "MJD: " << ev.mjd << " RA: " << ev.ra << " Dec: " << ev.dec << '\n';
    stream << "RA err: " << ev.ra_err_high << " " << ev.ra_err_low << '\n';
    stream << "Dec err: " << ev.dec_err_high << " " << ev.dec_err_low << '\n';

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const IceCubeEvent2& ev)
{
    stream << "ID: " << ev.id << '\n';
    stream << "MJD: " << ev.mjd << " RA: " << ev.ra << " Dec: " << ev.dec << '\n';
    stream << "Resolution: " << ev.resolution << " LogEnergy: " << ev.logEnergy << '\n';

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const AntaresEvent& ev)
{
    stream << "ID: " << ev.id << '\n';
    stream << "MJD: " << ev.mjd << " RA: " << ev.ra << " Dec: " << ev.dec << '\n';
    stream << "Nhits: " << ev.nhits << " Resolution: " << ev.resolution << '\n';

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const BaikalEvent& ev)
{
    stream << "ID: " << ev.id << '\n';
    stream << "MJD: " << ev.mjd << " RA: " << ev.ra << " Dec: " << ev.dec << '\n';
    stream << "Nhits: " << ev.nhits << " Resolution: " << ev.resolution << '\n';
    stream << "Energy: " << ev.energy << " EnergySigma: " << ev.energySigma << '\n';

    return stream;
}

std::ostream& operator<<(std::ostream& stream, const Coincidence& c)
{
    stream << "Events: " << c.det[0] << c.id << ", " << c.det[1] << c.id2 << '\n';
    stream << "timeDev: " << c.timeDev*24 << " angleDev: " << c.angleDev << '\n';

    return stream;
}

//binary search function for finding index of first event that happened later than
//given modified julian date (assuming events sorted by time)
int firstLaterIndex(const std::vector<Event>& arr, const double& mjd, const int& l = 0)
{
    //boundaries for search, current element is in the middle
    int low  = l;
    int high = arr.size();
    int middle;

    //searching index
    while (low <= high)
    {
        middle = (low + high) / 2;

        if (arr[middle].mjd < mjd)  low = middle + 1;
        else                       high = middle - 1;
    }
    
    if (arr[middle].mjd < mjd) middle++;

    if (middle >= arr.size()) middle = -1; //return -1 if their is no later event

    return middle;
}

//function for splitting string by whitespace into substrings 
std::vector<std::string> split(const std::string& str, char separator = ' ')
{
    std::vector<std::string> output;
    int currentPos = 0;
    int lastSpacePos = -1;

    for (char c : str)
    {
        if (c == separator) 
        {
            std::string substr = str.substr(lastSpacePos + 1, currentPos - lastSpacePos - 1);
            lastSpacePos = currentPos;
            if (substr != "") output.push_back(substr);
        }

        currentPos++;
    }

    std::string substr = str.substr(lastSpacePos + 1, currentPos - lastSpacePos - 1);
    if (substr != "") output.push_back(substr);

    return output;
}

//loading IceCube data from txt catalog of alerts
std::vector<IceCubeEvent> readIceCube(const char* path, int startID = 0)
{
    std::vector<std::string> textIn; //input text stripped of full line comments
    std::vector<IceCubeEvent> dataOut;
    std::ifstream inf{ path };
    int id = startID;

    //loading text file into dynamic string array
    while (inf)
    {
        std::string line;
        std::getline(inf, line);

        //ignoring comments, at(0) throws exception for empty string
        if(line != "" and line.at(0) != '#')
        {
            textIn.push_back(line);
        }
    }

    std::vector<std::string> textIn2; //input text striped of all comments

    for (std::string line : textIn) //ignoring comments on the end of line
    {
        int pos = line.find_first_of('#');
        if (pos != -1) textIn2.push_back(line.substr(0, pos));
        else textIn2.push_back(line);
    }

    //parsing each string into an IceCubeEvent object
    for (std::string line : textIn2)
    {
        std::vector<std::string> data = split(line);

        //saving each number from text into corresponding variable
        double mjd = std::stod(data[0]);
        double ra  = std::stod(data[1]);
        double dec = std::stod(data[3]);
        double ra_err_high = 0, ra_err_low = 0, dec_err_high = 0, dec_err_low = 0;

        //parsing errors if supplied
        if(data[2] != "(-)")
        {
            int commaIndex  = data[2].find_first_of(',');
            int commaIndex2 = data[4].find_first_of(',');

            int size  = data[2].size();
            int size2 = data[4].size();

            ra_err_high  = stod(data[2].substr(2, commaIndex  - 2));
            dec_err_high = stod(data[4].substr(2, commaIndex2 - 2));

            ra_err_low  = stod(data[2].substr(commaIndex  + 1, size  - commaIndex  - 2));
            dec_err_low = stod(data[4].substr(commaIndex2 + 1, size2 - commaIndex2 - 2));
        }

        dataOut.push_back(IceCubeEvent( 'I', id, mjd, ra, dec, ra_err_high, ra_err_low, dec_err_high, dec_err_low ));
        id++;
    }

    return dataOut;
}

//loading IceCube data from txt files (different formating)
std::vector<IceCubeEvent2> readIceCube2(const char* path, int startID = 0)
{
    std::vector<IceCubeEvent2> dataOut;
    std::ifstream inf{ path };
    int id = startID;

    //skipping first line
    {
        std::string line;
        std::getline(inf, line);
    }

    //reading and parsing entire file
    while (inf)
    {
        std::string Mjd, Ra, Dec, Resolution, LogEnergy;

        inf >> Mjd;
        inf >> Ra;
        inf >> Dec;
        inf >> Resolution;
        inf >> LogEnergy;

        if(Mjd != "")
        {
            double mjd          = stod(Mjd);
            double ra           = stod(Ra);
            double dec          = stod(Dec);
            double resolution   = stod(Resolution);
            double logEnergy    = stod(LogEnergy);

            dataOut.push_back(IceCubeEvent2( 'I', id, mjd, ra, dec, resolution, logEnergy ));
            id++;
        }
    }

    return dataOut;    
}

//loading IceCube data from new csv catalog of alerts
std::vector<IceCubeEvent> readIceCube3(const char* path, int startID = 0)
{
    std::vector<IceCubeEvent> dataOut;
    std::ifstream inf{ path };
    int id = startID;

    //skipping first line
    {
        std::string line;
        std::getline(inf, line);
    }

    //loading text file into dynamic string array
    while (inf)
    {
        std::string line;
        std::getline(inf, line);

        if(line != "")
        {
            std::vector<std::string> data = split(line, ',');

            if(data[1] != "NaN")
            {
                double ra  = std::stod(data[1]);
                double dec = std::stod(data[2]);

                double ra_err_high = 0, ra_err_low = 0, dec_err_high = 0, dec_err_low = 0;

                //parsing asymmetric errors
                if(data[3] != "NaN")
                {
                    int spaceIndex  = data[3].find_first_of(' ');
                    int spaceIndex2 = data[4].find_first_of(' ');

                    int size  = data[3].size();
                    int size2 = data[4].size();

                    ra_err_high  = std::stod(data[3].substr(1, spaceIndex  - 1));
                    dec_err_high = std::stod(data[4].substr(1, spaceIndex2 - 1));

                    ra_err_low  = std::stod(data[3].substr(spaceIndex  + 1, size  - spaceIndex  - 1));
                    dec_err_low = std::stod(data[4].substr(spaceIndex2 + 1, size2 - spaceIndex2 - 1));
                }

                //parsing time
                int year       = std::stoi(data[5].substr(0,4));
                int month      = std::stoi(data[5].substr(5,2));
                int day        = std::stoi(data[5].substr(8,2));
                int hours      = std::stoi(data[5].substr(11,2));
                int minutes    = std::stoi(data[5].substr(14,2));
                double seconds = std::stod(data[5].substr(17,2));

                double mjd = UTCtoMJD(year, month, day, hours, minutes, seconds);

                dataOut.push_back(IceCubeEvent('I',id,mjd,ra,dec,ra_err_high,ra_err_low,
                                               dec_err_high,dec_err_low));
                id++;
            }
        }
    }

    std::reverse(dataOut.begin(),dataOut.end());

    return dataOut;
}

std::vector<AntaresEvent> readAntares(const char* path, int startID = 0)
{
    std::vector<AntaresEvent> dataOut;
    std::ifstream inf{ path };
    int id = startID;

    //skipping first two lines
    {
        std::string line;
        std::getline(inf, line);
        std::getline(inf, line);
    }

    //reading and parsing entire file
    while (inf)
    {
        std::string Dec, Ra, Nhits, Resolution, Mjd;

        inf >> Dec;
        inf >> Ra;
        inf >> Nhits;
        inf >> Resolution;
        inf >> Mjd;

        if(Dec != "")
        {
            double dec          = stod(Dec);
            double ra           = stod(Ra);
            int nhits           = stoi(Nhits);
            double resolution   = stod(Resolution);
            double mjd          = stod(Mjd);

            dataOut.push_back(AntaresEvent( 'A', id, mjd, ra, dec, nhits, resolution ));
            id++;
        }
    }

    return dataOut;
}

std::vector<BaikalEvent> readBaikal(const char* path, int startID = 0)
{
    std::vector<BaikalEvent> dataOut;

    //access tree with data inside root file
    TFile data(path);
    TTree* cascades;

    data.GetObject("filteredCascades", cascades);
    
    //variables for data
    int id = startID;
    double ra, dec, resolution;
    int nhits;
    double e, eS;
    TTimeStamp* time = new TTimeStamp();

    
    cascades->SetBranchAddress("eventTime",       &time);
    cascades->SetBranchAddress("rightAscension",  &ra);
    cascades->SetBranchAddress("declination",     &dec);
    cascades->SetBranchAddress("directionSigma",  &resolution);
    cascades->SetBranchAddress("nHits",           &nhits);
    cascades->SetBranchAddress("energy",          &e);
    cascades->SetBranchAddress("energySigma",     &eS);

    //parsing data to BaikalEvent objects (converting time to MJD, angles to degrees)
    for (int i = 0; i < cascades->GetEntries(); ++i)
    {
        cascades->GetEntry(i);
        dataOut.push_back(BaikalEvent('B',id,time->AsJulianDate()-2400000.5,
                          radToDeg(ra),radToDeg(dec),resolution,nhits,e,eS));
        id++;
    }

    delete time;

    return dataOut;
}

//function for finding coincidences between events recorded on different experiments
//input = 2 vectors of events (obtained by read functions), maximal time difference (days),
//maximal angle difference (degrees), time shift of first set (days)
//second vector has to be sorted by time of the event (mjd variable)

//passing smaller vector as first argument may increase the speed significantly
std::vector<Coincidence> findCoincidences(const std::vector<Event>& events, const std::vector<Event>& events2,
                                          const double& maxTimeDev, const double& maxAngleDev,
                                          const double& timeShift = 0)
{
    std::vector<Coincidence> output;

    std::string det; det += events[0].det; det += events2[0].det; //detector codes
    
    for(auto ev : events)
    {
        double minTime = ev.mjd - maxTimeDev + timeShift;
        double maxTime = ev.mjd + maxTimeDev + timeShift;

        int i = firstLaterIndex({events2.begin(), events2.end()}, minTime);

        if(i != -1)
        { 
            for( ; events2[i].mjd < maxTime && i < events2.size(); i++)
            {
                double angleDev = angularDistance(ev.ra, ev.dec, events2[i].ra, events2[i].dec);

                if(angleDev < maxAngleDev)
                {
                    double timeDev = std::abs(ev.mjd - events2[i].mjd);
                    output.push_back(Coincidence{det, ev.id, events2[i].id, timeDev, angleDev});
                }
            }
        }
    }
    
    return output;
}

//checks if splitted array is ascending using first 4 elements
bool isAscending(const double* arr)
{
    int num = (arr[0] < arr[1]) + (arr[1] < arr[2]) + (arr[2] < arr[3]);
    return num > 1; 
}

//sort x and y array with splitted data
template<int size>
XYarr<size> sortSplitted(const double* x, const double* y)
{
    int splitIndex = 0;

    if(isAscending(x))
    {
        for (int i = 0; i < size-1; ++i) 
            if(x[i] > 0 and x[i+1] < 0) splitIndex = i+1;
    }

    else
    {
        for (int i = 0; i < size-1; ++i) 
            if(x[i] < 0 and x[i+1] > 0) splitIndex = i+1;       
    }

    double x2[size], y2[size];

    for (int i = splitIndex; i < size; ++i)
    {
        x2[i - splitIndex] = x[i];
        y2[i - splitIndex] = y[i];
    }

    for (int i = 0; i < splitIndex; ++i)
    {
        x2[i + size - splitIndex] = x[i];
        y2[i + size - splitIndex] = y[i];
    }

    return XYarr<size>{x2, y2};
}

//function for plotting random coincidences into histogram
//banned window width describes band around zero time shift where real coincidences might occur
//maxShift on both sides, banWinWidth is band width
TH1F* randomCoincidences(const std::vector<Event>& events, const std::vector<Event>& events2,
                         const double& maxTimeDev, const double& maxAngleDev, const int& samples,
                         const double& maxShift, const double& banWinWidth)
{
    TRandom2 random(0);

    double numOfCoincidences[samples];

    for(int i = 0; i < samples; ++i)
    {
        double shift;

        do shift = 2*maxShift*(random.Rndm()-0.5);
        while(std::abs(shift) < banWinWidth/2);

        std::vector<Coincidence> c = findCoincidences(events, events2, maxTimeDev, maxAngleDev, shift);

        numOfCoincidences[i] = c.size();
    }

    double maxC = maxElement(numOfCoincidences, samples);
    TH1F* histogram = new TH1F("histogram", "Random coincidences", maxC, 0, maxC);

    histogram->FillN(samples, numOfCoincidences, nullptr);
    histogram->Scale(1./samples);
    histogram->Draw("hist");

    return histogram;
}

int coincidence()
{
    std::vector<IceCubeEvent> iceEv = readIceCube("experiments_data/catalog_of_alerts.txt");
    std::vector<AntaresEvent> antEv = readAntares("experiments_data/antares_events.txt");
    std::vector<BaikalEvent>  baiEv = readBaikal("experiments_data/filteredCascades_y18c-1.root");

    std::vector<IceCubeEvent2> ice2Ev = readIceCube2("experiments_data/IceCube_data/events_IC40.txt");

    ice2Ev += readIceCube2("experiments_data/IceCube_data/events_IC59.txt",  ice2Ev.size());
    ice2Ev += readIceCube2("experiments_data/IceCube_data/events_IC79.txt",  ice2Ev.size());
    ice2Ev += readIceCube2("experiments_data/IceCube_data/events_IC86a.txt", ice2Ev.size());
    ice2Ev += readIceCube2("experiments_data/IceCube_data/events_IC86b.txt", ice2Ev.size());
    ice2Ev += readIceCube2("experiments_data/IceCube_data/events_IC86c.txt", ice2Ev.size());

    std::vector<IceCubeEvent> iceEv2 = readIceCube3("experiments_data/events.csv");

    randomCoincidences({baiEv.begin(), baiEv.end()}, {iceEv2.begin(),iceEv2.end()}, 1, 20, 10000, 50, 2);

    std::vector<Coincidence> coincidences = findCoincidences({baiEv.begin(), baiEv.end()}, {iceEv2.begin(),iceEv2.end()}, 1, 20);

    for(auto ev : coincidences) std::cout << ev << '\n';

    //std::cout << ice2Ev.size() << '\n';

    //for(auto ev : iceEv2) std::cout << ev << '\n';

    /*TCanvas* c1 = new TCanvas("c1", "", 1000, 500);
    drawmap("Coincidences between Baikal-GVD and ANTARES;Right ascension;Declination");

    TGraph* gBai = new TGraph();
    TGraph* gAnt = new TGraph();

    for(Coincidence c : coincidences)
    {
        double rAsc = baiEv[c.id].ra;
        if(rAsc>180) rAsc -= 360;

        double rAsc2 = antEv[c.id2].ra;
        if(rAsc2>180) rAsc2 -= 360;

        XY baiCoor = toAitoff(rAsc, baiEv[c.id].dec);
        XY antCoor = toAitoff(rAsc2, antEv[c.id2].dec);

        gBai->SetPoint(gBai->GetN(),baiCoor.x, baiCoor.y);
        gAnt->SetPoint(gAnt->GetN(),antCoor.x, antCoor.y);
    }

    gBai->SetMarkerStyle(kCircle);
    gAnt->SetMarkerStyle(kStar);

    gBai->SetMarkerColor(kBlue);
    gAnt->SetMarkerColor(kRed);

    gBai->Draw("P");
    gAnt->Draw("P");

    drawLabels();*/
    
    return 0;
}

    /*for(AntaresEvent ev : antEv)
    {
        double rAsc = ev.ra;
        if(rAsc>180) rAsc -= 360;
        XY baiCoor = toAitoff(rAsc, ev.dec);

        gBai->SetPoint(gBai->GetN(),baiCoor.x, baiCoor.y);
    }*/