/* 
 * File	  : Timer.hpp 
 * Author : Md Hasanuzzaman Bhuiyan
 * Email  : mhb@vbi.vt.edu
 * Date	  : July 27, 2012
 */

#ifndef _Timer_HPP
#define	_Timer_HPP

#include <ctime>
#include <sys/time.h>

class Timer
{
	struct timeval starttime, endtime, diff;
	long diff_microsec;	// difference of time is saved in this variable in microsecond format
	bool s;		// s = true means start function was called, false means start function was not started
	bool e;		// e = true means end function was called, false means end function was not called

	/*------- Calculate the difference of two time in "micro-second" ----------*/	
	void TimeDifference()
	{
		diff.tv_sec  = endtime.tv_sec  - starttime.tv_sec;
		diff.tv_usec = endtime.tv_usec - starttime.tv_usec;
		while (diff.tv_usec < 0)
		{
			diff.tv_usec += 1000000;
			diff.tv_sec -= 1;
		}
			
		diff_microsec = 1000000L * diff.tv_sec + diff.tv_usec;	// returning time difference in "micro-second"
	}
	
	public:
		
		Timer() 	// constructor : do nothing
		{ 
			s = e = false; 
		}
		
		Timer(int i)	// constructor taking an integer parameter : starting time count
		{  
			gettimeofday(&starttime, NULL);
			s = true;
			e = false;
		}
		
		void start ()	// start function : starting time count (called when constructor was called without parameter)
		{
			gettimeofday(&starttime, NULL);
			s = true;
			e = false;
		}
		
		// ends time counting and return the time difference in microseconds (return data type: long)
		long end()
		{
			if(s)
			{
				gettimeofday(&endtime, NULL);
				TimeDifference();
				e = true;
				return diff_microsec;
			}
			else
			{
				printf("Error!\n\"end\" function is called without calling the \"start\" function\n");
				return -1;
			}
		}
		
		// prints time difference in micro second
		void printmicsec()
		{
			if(!s)
			{
				printf("Error!\n\"print\" function is called without calling the \"start\" function\n");
				return;
			}
			if(!e)
				end();
			printf("%ld micro-second", diff_microsec);			
		}
		
		// prints time difference in milli second
		void printmilsec()
		{
			if(!s)
			{
				printf("Error!\n\"print\" function is called without calling the \"start\" function\n");
				return;
			}
			if(!e)
				end();
			printf("%.2lf milli-second", ((double)diff_microsec)/1000.0);			
		}
		
		// prints time difference in second
		void printsec()
		{
			if(!s)
			{
				printf("Error!\n\"print\" function is called without calling the \"start\" function\n");
				return;
			}
			if(!e)
				end();
			printf("%.2lf second", ((double)diff_microsec)/1000000.0);			
		}
		
		double getsec()
		{
			if(!s)
			{
				printf("Error!\n\"print\" function is called without calling the \"start\" function\n");
				return 0.0L;
			}
			if(!e)
				end();
			
			return ((double)diff_microsec)/1000000.0;
		}
		
		// prints time difference in minute and second
		void printminsec()
		{
			if(!s)
			{
				printf("Error!\n\"print\" function is called without calling the \"start\" function\n");
				return;
			}
			if(!e)
				end();
			long minute = diff_microsec/60000000;
			printf("%ld minute %.2lf second", minute, (diff_microsec/1000000.0) - minute*60.0);			
		}
		
		// prints time difference in hour minute second
		void printhourminsec()
		{
			if(!s)
			{
				printf("Error!\n\"print\" function is called without calling the \"start\" function\n");
				return;
			}
			if(!e)
				end();
			long hour = (diff_microsec/60000000)/60;
			long minute = (diff_microsec/60000000)%60;
			printf("%ld hour %ld minute %.2lf second", hour, minute, ((double)diff_microsec)/1000000.0 - hour*3600.0 - minute*60.0);			
		}
};

#endif /* _Timer_HPP */
