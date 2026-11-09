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
    TFile file(fileName, "UPDATE");
    if (file.IsZombie()) {
        printf("Failed to open %s\n", fileName);
        return;
    }

    std::set<std::string> directoryNames;
    TIter nextKey(file.GetListOfKeys());

    while (auto* key = dynamic_cast<TKey*>(nextKey())) {
        TString name = key->GetName();

        if (name.BeginsWith("multiparticle-cumulants")) {
            directoryNames.insert(name.Data());
        }
    }

    Long64_t totalRemoved = 0;
    Long64_t totalKept = 0;

    for (const std::string& directoryName : directoryNames) {
        TDirectory* directory =
            file.GetDirectory(directoryName.c_str());

        if (!directory) {
            continue;
        }

        TList* baseList = nullptr;
        directory->GetObject("Default list name", baseList);

        if (!baseList) {
            printf("%s: Default list not found\n",
                   directoryName.c_str());
            continue;
        }

        CleanCounts counts;
        CleanCollection(baseList, removeAll, counts);

        directory->cd();

        const int bytesWritten = baseList->Write(
            "Default list name",
            TObject::kSingleKey | TObject::kWriteDelete);

        if (bytesWritten <= 0) {
            printf("%s: write failed\n",
                   directoryName.c_str());
            continue;
        }

        totalRemoved += counts.removed;
        totalKept += counts.kept;

        printf("%s: removed %lld TF1, kept %lld\n",
               directoryName.c_str(),
               counts.removed,
               counts.kept);
    }

    file.Close();

    printf("Total removed: %lld TF1, kept: %lld\n",
           totalRemoved,
           totalKept);
}