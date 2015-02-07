#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED


#include <vector>


#ifdef WIN32
	#include <windows.h>
	#define time_type LARGE_INTEGER
	#define MAC_getTime(X) QueryPerformanceCounter(&X);
	
	#define MAC_addDifference(dest, a, b) dest.QuadPart += a.QuadPart - b.QuadPart;
#else
	#include <sys/time.h>
	#define time_type timeval
	#define MAC_getTime(X) gettimeofday(&X,0)

	#define MAC_addDifference(dest, a, b) dest.tv_sec += a.tv_sec - b.tv_sec; dest.tv_usec += a.tv_usec - b.tv_usec;
#endif //WIN32


namespace {
	
	double toSeconds(const time_type& totalTime) {
#ifdef WIN32
		LARGE_INTEGER freq; 
		QueryPerformanceFrequency(&freq);
		return (double)totalTime.QuadPart / (double)freq.QuadPart;
#else
		return (double)totalTime.tv_sec + totalTime.tv_usec*0.000001;
#endif //WIN32
	}
	
	double toMicroSeconds(const time_type& totalTime) {
#ifdef WIN32
		LARGE_INTEGER freq; 
		QueryPerformanceFrequency(&freq);
		return (double)totalTime.QuadPart * 1000000.0 / (double)freq.QuadPart;
#else
		return (double)totalTime.tv_sec*1000000.0 + totalTime.tv_usec;
#endif //WIN32
	}
	
}


class PrecisionTimer1 {
	time_type startTime;
	time_type totalTime;

public:
	PrecisionTimer1() {
		time_type temp = {0};
		totalTime = temp;
	}
	
	void start() {
		MAC_getTime(startTime);
	}
	
	time_type read() {
		time_type temp;
		MAC_getTime(temp);
		return temp;
	}
	
	void pause() {
		time_type endTime;
		MAC_getTime(endTime);
		MAC_addDifference(totalTime, endTime, startTime);
	}

	double totalSeconds() const {
		return toSeconds(totalTime);
	}
	
	operator time_type() {
		return totalTime;
	}
};


class PrecisionTimerSequence {
	std::vector<std::pair<time_type, int> > sequence;
	
public:
	void reserve(size_t rsize) {
		sequence.reserve(rsize);
	}
	
	void addTime(const time_type& t, int id) {
		sequence.push_back(std::make_pair(t, id));
	}
	
	size_t size() const {
		return sequence.size();
	}
	
	std::pair<double, int> at(size_t i) const {
		time_type temp = {0};
		MAC_addDifference(temp, sequence[i].first, sequence[0].first);
		return std::make_pair(toSeconds(temp), sequence[i].second);
	}
};


class PrecisionTimer {
	time_type startTime;
	time_type totalTime;
	PrecisionTimerSequence* seq;
	int id;
	bool paused;

public:
	PrecisionTimer() : seq(0), paused(true) {
		time_type temp = {0};
		totalTime = temp;
	}
	
	PrecisionTimer(PrecisionTimerSequence& s, int i) : seq(&s), id(i), paused(true) {
		time_type temp = {0};
		totalTime = temp;
	}
	
	void start() {
	    if (paused) {
            paused = false;
            MAC_getTime(startTime);
            if (seq)
                seq->addTime(startTime, id);
	    } else
            throw "not paused when trying to start";
	}
	
	time_type read() {
		time_type temp;
		MAC_getTime(temp);
		return temp;
	}
	
	void pause() {
	    if (!paused) {
	        paused = true;
            time_type endTime;
            MAC_getTime(endTime);
            MAC_addDifference(totalTime, endTime, startTime);
            if (seq)
                seq->addTime(endTime, -id);
	    } else throw "already paused when trying to pause";
	}

	double totalSeconds() const {
	    if (!paused)
            throw "timer not paused when total number of seconds is being read!";
		return toSeconds(totalTime);
	}
	
	operator time_type() {
		return totalTime;
	}
	
	void clear() {
	    paused = true;
		time_type temp = {0};
		totalTime = temp;
	}
	
	/// wait until total time passed is larger than or equal to "tm"
	void waitTotal(double tm) {
	    //do {MAC_getTime(t); MAC_addDifference(tt, t, startTime);} while (toMicroSeconds(tt) < tm);
	    do {
	        time_type t, tt = {0};
	        MAC_getTime(t);
	        MAC_addDifference(tt, t, startTime);
	        if (toSeconds(tt) >= tm*0.000001)
                break;
	    } while(true);
	}
};


// create a scope timer to add the duration of current scope to selected PrecisionTimer
class ScopeTimer {
	PrecisionTimer& pt;
	
public:
	ScopeTimer(PrecisionTimer& t) : pt(t) {
		if (&pt != 0) pt.start();
	}
	
	~ScopeTimer() {
		if (&pt != 0) pt.pause();
	}
};


class ExcludeScopeFromTimer {
	PrecisionTimer& pt;
	
public:
	ExcludeScopeFromTimer(PrecisionTimer& t) : pt(t) {
		if (&pt != 0) pt.pause();
	}
	
	~ExcludeScopeFromTimer() {
		if (&pt != 0) pt.start();
	}
};



#endif // TIMER_H_INCLUDED
