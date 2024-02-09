#include "refcontribwriter.h"

namespace raysar
{

// static instance member needs to be declared outside the class body to avoid LNK2001
boost::shared_ptr<RefContribWriter> RefContribWriter::mInstance;

boost::shared_ptr<RefContribWriter> RefContribWriter::getInstance()
{
   // lock_guard is not very efficient here, but since we can't expect C++11 features, there are no Magic Statics available
   boost::lock_guard<boost::mutex> myGuard(mMutex);
   if (!mInstance)
   {
      // Still no singleton instance available, an instance needs to be created
      mInstance = boost::make_shared<RefContribWriter>(mFilename);
   }

   return mInstance;
}

void RefContribWriter::add(DBL azimuth, DBL slantrange, DBL elevation, DBL intensity, 
            int tracelevel, int specular, DBL sx, DBL sy, DBL sz)
{
   auto contrib = boost::make_shared<RefContrib>(azimuth, slantrange, elevation,
                                 intensity, tracelevel, specular, sx, sy, sz);

   {
      boost::lock_guard<boost::mutex> myGuard(mMutex);
      // Add contribution to the contributions list
      mContribs.push_front(std::move(contrib));
      mIntersectFlags.push_front(1);
   }
}

void RefContribWriter::add(DBL azimuth, DBL slantrange, DBL elevation, DBL intensity, 
            int tracelevel, int specular)
{
   auto contrib = boost::make_shared<RefContrib>(azimuth, slantrange, elevation,
                                 intensity, tracelevel, specular, 0.0, 0.0, 0.0);

   {
      boost::lock_guard<boost::mutex> myGuard(mMutex);
      // Add contribution to the contributions list
      mContribs.push_front(std::move(contrib));
      mIntersectFlags.push_front(0);
   }
}

void RefContribWriter::finish()
{
   stop();
   mWriterThread.join();
}

void RefContribWriter::writeToFile()
{
   boost::lock_guard<boost::mutex> myGuard(mMutex);
   
   while (!mContribs.empty())
   {
      boost::shared_ptr<RefContrib> myContrib = mContribs.back();
      char myIntersectFlag = mIntersectFlags.back();
      mContribs.pop_back();
      mIntersectFlags.pop_back();
      
      if (myIntersectFlag > 0)
      {
         fprintf(mFid, "%.6f %.6f %.6f %.6f %d %d %.6f %.6f %.6f \n",  
                  myContrib->mAzimuth,
                  myContrib->mSlantRange,
                  myContrib->mElevation,
                  myContrib->mIntensity,
                  myContrib->mTraceLevel,
                  myContrib->mSpecular,
                  myContrib->mSX,
                  myContrib->mSY,
                  myContrib->mSZ);
      }
      else 
      {
         fprintf(mFid, "%.6f %.6f %.6f %.6f %d %d \n",  
                  myContrib->mAzimuth,
                  myContrib->mSlantRange,
                  myContrib->mElevation,
                  myContrib->mIntensity,
                  myContrib->mTraceLevel,
                  myContrib->mSpecular);
      }
   }
}

void RefContribWriter::writeLoop()
{
   // setlocale(LC_ALL, "");
   setlocale(LC_NUMERIC, "en-US");

   // Open file
   mFid = fopen(mFilename, "w");
#ifdef RCW_DEBUG
   printf("RefContribWriter, %s opened, mFid created.\n", mFilename);
#endif

   while(isRunning())
   {
        // Continuously write contents to file
        writeToFile();
   }
   
   // Clear contribs and intersect save flags
   // (they should be empty anyhow)
   mContribs.clear();
   mIntersectFlags.clear();
   
   // Close file
   fclose(mFid);
#ifdef RCW_DEBUG
   printf("RefContribWriter, mFid closed.\n");
#endif
}

RefContribWriter::~RefContribWriter()
{
   finish();
#ifdef RCW_DEBUG
   printf("RefContribWriter::~RefContribWriter()\n");
#endif
}

}
