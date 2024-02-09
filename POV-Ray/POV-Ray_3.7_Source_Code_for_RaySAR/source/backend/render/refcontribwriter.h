#ifndef _REFCONTRIBS_H_
#define _REFCONTRIBS_H_

#include <deque>
#include <memory>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/make_shared.hpp>
#include <algorithm>
#include <cstdio>
#include <locale>

#ifndef DBL
   typedef double DBL;
#endif

// #define RCW_DEBUG 1

namespace raysar 
{
   struct RefContrib
   {
   public:
      DBL mAzimuth;
      DBL mSlantRange;
      DBL mElevation;
      DBL mIntensity; 
      int mTraceLevel;
      int mSpecular;
      DBL mSX;
      DBL mSY;
      DBL mSZ;

      RefContrib(DBL azimuth, DBL slantrange, DBL elevation, DBL intensity, 
               int tracelevel, int specular, DBL sx, DBL sy, DBL sz) :
         mAzimuth(azimuth),
         mSlantRange(slantrange),
         mElevation(elevation),
         mIntensity(intensity), 
         mTraceLevel(tracelevel),
         mSpecular(specular),
         mSX(sx),
         mSY(sy),
         mSZ(sz)
      {
      }
   };

   class RefContribWriter
   {
   protected:
      static boost::shared_ptr<RefContribWriter>   mInstance;

      std::deque<boost::shared_ptr<RefContrib> >   mContribs;
      std::deque<char>                             mIntersectFlags;
      boost::mutex                                 mMutex;
      boost::mutex                                 mRunMutex;

      boost::thread                                mWriterThread;
      bool                                         mRunning;
      const char*                                  mFilename;
      std::FILE*                                   mFid;

   public:
      RefContribWriter(const char* filename) :
         mFilename(filename)
      {
#ifdef RCW_DEBUG
         printf("RefContribWriter, leaving constructor.\n");
#endif
      }

      // Disable copy construction for singleton
      RefContribWriter(RefContribWriter const& r);		// Do not implement
      RefContribWriter& operator=(RefContribWriter const& r); 	// Do not implement

      ~RefContribWriter();

      boost::shared_ptr<RefContribWriter> getInstance();
      
      void add(DBL azimuth, DBL slantrange, DBL elevation, DBL intensity, 
               int tracelevel, int specular, DBL sx, DBL sy, DBL sz);
      void add(DBL azimuth, DBL slantrange, DBL elevation, DBL intensity, 
               int tracelevel, int specular);

      void finish();

      // Basic start / stop
      void start()
      {
         boost::lock_guard<boost::mutex> myGuard(mRunMutex);
         mRunning = true;
#ifdef RCW_DEBUG
         printf("RefContribWriter, creating thread ...\n");
#endif
         boost::thread t (boost::bind(&RefContribWriter::writeLoop, this));
         mWriterThread.swap(t);
      }
      
      void stop()
      {
         boost::lock_guard<boost::mutex> myGuard(mRunMutex);
         mRunning = false;
      }

   protected:
      void writeToFile();
      void writeLoop();
      
      bool isRunning()
      {
         bool res;
         
         {
            boost::lock_guard<boost::mutex> myGuard(mRunMutex);
            res = mRunning;
         }
         
         return res;
      }
      
   };

}

#endif
