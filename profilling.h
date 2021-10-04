#pragma once

#if PROFILLING
	#include <chrono>
	#include <string>
	#include <vector>

	class AdditiveTimer
	{
	private:
		const char* m_Name;
		double m_Time = 0.0;
	friend class Profiler;

	public:
		AdditiveTimer(const char* name) : m_Name(name) {}
		//~AdditiveTimer();
		
		void AddToTimer(double newTime) {m_Time += newTime;}
	};

	class Profiler
	{
	private:
		std::vector<AdditiveTimer> m_Timers;

	public:
		Profiler(const Profiler&) = delete;

		AdditiveTimer* GetTimer(const char* name)
		{
			for(int i = 0; i < m_Timers.size(); i++) 
			{	
				if(std::strcmp(m_Timers[i].m_Name,name) == 0) return &m_Timers[i];
			}

			m_Timers.emplace_back(name);

			return &m_Timers.back();
		}

		double LastTime()
		{
			return m_Timers.back().m_Time;
		}

		static Profiler& Get()
		{
			static Profiler instance;
			return instance;
		}

	private:
		Profiler(){}
		//~Profiler();
	};

	//timer times either the entire scope or until Stop() is called
	class ProfillingTimer
	{
	private:
	    std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
	    const char* m_Name;
	    bool m_Stopped;

	public:
	    ProfillingTimer(const char* name) : m_Name(name), m_Stopped(false)
	    {
	        m_StartTimepoint = std::chrono::high_resolution_clock::now();
	    }

	    ~ProfillingTimer()
	    {
	        if (!m_Stopped) Stop();
	    }

	    void Stop()
	    {
	        auto endTimepoint = std::chrono::high_resolution_clock::now();

	        std::chrono::duration<double> duration = endTimepoint - m_StartTimepoint;
	        m_Stopped = true;

	        //std::cout << "Timer " << m_Name << " took " << duration.count() << "s\n";

	        Profiler::Get().GetTimer(m_Name)->AddToTimer(duration.count());
	    }
	};

#endif //PROFILLING