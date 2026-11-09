#include <TCollection.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TKey.h>
#include <TList.h>
#include <TString.h>

#include <set>
#include <string>
#include <vector>
#include <TROOT.h>
#include <algorithm>
#include <mutex>
#include <thread>

struct CleanCounts {
    Long64_t removed = 0;
    Long64_t kept = 0;
};

void CleanObject(TObject* object, bool removeAll, CleanCounts& counts);

void CleanCollection(TCollection* collection,
                     bool removeAll,
                     CleanCounts& counts)
{
    if (!collection) {
        return;
    }

    std::set<std::string> keptNames;
    std::vector<TObject*> objectsToRemove;

    TIter next(collection);
    while (TObject* object = next()) {
        if (object->InheritsFrom(TF1::Class())) {
            const bool duplicate =
                !keptNames.insert(object->GetName()).second;

            if (removeAll || duplicate) {
                objectsToRemove.push_back(object);
            } else {
                ++counts.kept;
            }
        } else {
            // 處理內層的 TList/TCollection 或 histogram
            CleanObject(object, removeAll, counts);
        }
    }

    for (TObject* object : objectsToRemove) {
        collection->Remove(object);
        delete object;
        ++counts.removed;
    }
}

void CleanObject(TObject* object,
                 bool removeAll,
                 CleanCounts& counts)
{
    if (!object) {
        return;
    }

    // 巢狀的 TList、TObjArray 等
    if (auto* collection = dynamic_cast<TCollection*>(object)) {
        CleanCollection(collection, removeAll, counts);
        return;
    }

    // TH1F、TH2F、TProfile 等附帶的 TF1
    if (auto* histogram = dynamic_cast<TH1*>(object)) {
        CleanCollection(histogram->GetListOfFunctions(),
                        removeAll,
                        counts);
    }
}



void CleanTF1(const char* fileName, bool removeAll = false)
{

    ROOT::EnableThreadSafety();
    std::mutex filemutex;

    TFile file(fileName, "UPDATE");
    if (file.IsZombie()) {
        printf("Failed to open %s\n", fileName);
        return;
    }

    std::set<std::string> directoryNames_set;
    TIter nextKey(file.GetListOfKeys());

    while (auto* key = dynamic_cast<TKey*>(nextKey())) {
        TString name = key->GetName();

        if (name.BeginsWith("multiparticle-cumulants")) {
            directoryNames_set.insert(name.Data());
        }
    }

    Long64_t totalRemoved = 0;
    Long64_t totalKept = 0;

    std::vector<std::string> directoryNames(
        directoryNames_set.begin(),
        directoryNames_set.end()
    );
    size_t totalSize = directoryNames.size();
    size_t nThreads = 4;
    size_t threadSize = totalSize / nThreads;
    size_t threadSize_remainder = totalSize % nThreads;
    std::vector<std::thread> threads;


    auto LoopOverDirectory = [&](std::string directoryName, bool removeAll = false) {

        TDirectory* directory = nullptr;
        TList* baseList = nullptr;
        {   
            std::lock_guard<std::mutex> lock(filemutex);

            directory = file.GetDirectory(directoryName.c_str());

            if (!directory) {
                return 1;
            }
            directory->GetObject("Default list name", baseList);
        }

        if (!baseList) {
            printf("%s: Default list not found\n",
                    directoryName.c_str());
            return 2;
        }

        auto* correlationList = dynamic_cast<TCollection*>(
            baseList->FindObject("CorrelationHistograms"));

        if (!correlationList) {
            printf("%s: CorrelationHistograms not found\n",
                directoryName.c_str());
            return 3;
        }

        CleanCounts counts;

        CleanCollection(correlationList, removeAll, counts);

        directory->cd();

        {
            std::lock_guard<std::mutex> lock(filemutex);
            if (counts.removed == 0) {
                totalKept += counts.kept;
                printf("%s: no TF1 removed; skipped writing\n",
                    directoryName.c_str());
                return 4;
            }

            const int bytesWritten = baseList->Write(
                "Default list name",
                TObject::kSingleKey | TObject::kWriteDelete);

            if (bytesWritten <= 0) {
                printf("%s: write failed\n",
                        directoryName.c_str());
                return 5;
            }
        
            totalRemoved += counts.removed;
            totalKept += counts.kept;

            printf("%s: removed %lld TF1, kept %lld\n",
                    directoryName.c_str(),
                    counts.removed,
                    counts.kept);
        }

        return 0;

    };


    for (int thread_id = 0; thread_id < nThreads; thread_id++) {
        size_t local_threadSize = threadSize;
        if (thread_id == nThreads - 1) {
            local_threadSize += threadSize_remainder;
        }

        std::vector<std::string> directoryName;
        for (int dirname_id = thread_id * threadSize; dirname_id < thread_id * threadSize + local_threadSize; dirname_id++){
            directoryName.push_back(directoryNames[dirname_id]);
        }


        threads.emplace_back( [&, directoryName](){
            for (auto& name : directoryName) {
                LoopOverDirectory(name, removeAll);
            }
        });
    }
    
    for (int thread_id = 0; thread_id < nThreads; thread_id++) {
        threads[thread_id].join();
    }

    file.Close();

    printf("Total removed: %lld TF1, kept: %lld\n",
           totalRemoved,
           totalKept);
}