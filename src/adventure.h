#ifndef SRC_ADVENTURE_H_
#define SRC_ADVENTURE_H_

#include <algorithm>
#include <vector>

#include "../third_party/threadpool/threadpool.h"

#include "./types.h"
#include "./utils.h"

class Adventure {
 public:
  virtual ~Adventure() = default;

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag& bag) = 0;

  virtual void arrangeSand(std::vector<GrainOfSand>& grains) = 0;

  virtual Crystal selectBestCrystal(std::vector<Crystal>& crystals) = 0;

 protected:
  void merge(std::vector<GrainOfSand>& grains, int i1, int i2, int end) {
    std::vector<GrainOfSand> res;
    int midEnd = i2 - 1;
    int begin = i1;

    while (i1 <= midEnd || i2 <= end) {
      if (i1 > midEnd || (i2 <= end && grains[i2] < grains[i1])) {
        res.push_back(grains[i2]);
        i2++;
      } else {
        res.push_back(grains[i1]);
        i1++;
      }
    }

    for (int i = begin; i <= end; i++) grains[i] = res[i - begin];
  }

  void sort(std::vector<GrainOfSand>& grains, int left, int right) {
    if (left == right) return;

    int mid = (left + right) / 2;
    sort(grains, left, mid);
    sort(grains, mid + 1, right);
    merge(grains, left, mid + 1, right);
  }

  int fillTheBag(BottomlessBag& bag, std::vector<Egg>& eggs, int** dp, int cap,
                 int pos) {
    int size = eggs[pos].getSize();
    if (pos == 0 && dp[pos][cap] != 0) {
      bag.addEgg(eggs[0]);
    } else if (pos != 0) {
      if (dp[pos][cap] == dp[pos - 1][cap]) {
        fillTheBag(bag, eggs, dp, cap, pos - 1);
      } else {
        bag.addEgg(eggs[pos]);
        fillTheBag(bag, eggs, dp, cap - size, pos - 1);
      }
    }

    return dp[pos][cap];
  }
};

class LonesomeAdventure : public Adventure {
 public:
  LonesomeAdventure() {}

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag& bag) {
    int** dp = new int*[eggs.size()];
    for (unsigned int i = 0; i < eggs.size(); i++)
      dp[i] = new int[bag.getCapacity() + 1];

    for (unsigned int egg = 0; egg < eggs.size(); egg++) {
      for (unsigned int cap = 0; cap <= bag.getCapacity(); cap++) {
        int weight = eggs[egg].getWeight();
        unsigned int size = eggs[egg].getSize();

        if (egg == 0) {
          dp[egg][cap] = (size > cap) ? 0 : weight;
        } else {
          dp[egg][cap] = dp[egg - 1][cap];
          if (cap >= size)
            dp[egg][cap] =
                std::max(dp[egg][cap], dp[egg - 1][cap - size] + weight);
        }
      }
    }

    int res = fillTheBag(bag, eggs, dp, bag.getCapacity(), eggs.size() - 1);
    for (unsigned int i = 0; i < eggs.size(); i++) delete[] dp[i];
    delete[] dp;
    return res;
  }

  virtual void arrangeSand(std::vector<GrainOfSand>& grains) {
    sort(grains, 0, grains.size() - 1);
  }

  virtual Crystal selectBestCrystal(std::vector<Crystal>& crystals) {
    int maxPos = 0;
    for (unsigned int i = 0; i < crystals.size(); i++)
      if (crystals[maxPos] < crystals[i]) maxPos = i;

    return crystals[maxPos];
  }
};

class TeamAdventure : public Adventure {
 public:
  explicit TeamAdventure(uint64_t numberOfShamansArg)
      : numberOfShamans(numberOfShamansArg),
        councilOfShamans(numberOfShamansArg) {}

  uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag& bag) {
    int** dp = new int*[eggs.size()];
    for (unsigned int i = 0; i < eggs.size(); i++)
      dp[i] = new int[bag.getCapacity() + 1];

    auto intervals = getIntervals(bag.getCapacity(), numberOfShamans);
    std::vector<std::future<void>> futures;

    for (unsigned int i = 0; i < eggs.size(); i++) {
      futures.clear();
      for (auto interval : intervals)
        futures.push_back(
            councilOfShamans.enqueue([this, dp, &eggs, i, interval]() {
              knapsack(dp, eggs, i, interval.first, interval.second);
            }));

      for (auto& f : futures) f.wait();
    }

    int res = fillTheBag(bag, eggs, dp, bag.getCapacity(), eggs.size() - 1);
    for (unsigned int i = 0; i < eggs.size(); i++) delete[] dp[i];
    delete[] dp;
    return res;
  }

  virtual void arrangeSand(std::vector<GrainOfSand>& grains) {
    auto intervals = getIntervals(grains.size() - 1, numberOfShamans);
    std::vector<std::future<void>> futures;
    for (auto i : intervals)
      futures.push_back(councilOfShamans.enqueue(
          [this, &grains, i]() { sort(grains, i.first, i.second); }));

    for (auto& f : futures) f.wait();

    mergeIntervals(intervals, grains, 0, intervals.size() - 1);
  }

  virtual Crystal selectBestCrystal(std::vector<Crystal>& crystals) {
    int res = 0;
    std::mutex mut;
    auto intervals = getIntervals(crystals.size() - 1, numberOfShamans);
    std::vector<std::future<void>> futures;
    for (auto interval : intervals) {
      futures.push_back(
          councilOfShamans.enqueue([&crystals, &mut, &res, interval]() {
            int maxPos = interval.first;
            for (unsigned int i = interval.first; i <= interval.second; i++)
              if (crystals[maxPos] < crystals[i]) maxPos = i;

            mut.lock();
            if (crystals[res] < crystals[maxPos]) res = maxPos;
            mut.unlock();
          }));
    }

    for (auto& f : futures) f.wait();

    return crystals[res];
  }

 private:
  uint64_t numberOfShamans;
  ThreadPool councilOfShamans;
  using Interval = std::pair<unsigned int, unsigned int>;

  void knapsack(int** dp, std::vector<Egg>& eggs, int egg, int left,
                int right) {
    for (int cap = left; cap <= right; cap++) {
      int size = eggs[egg].getSize();
      int weight = eggs[egg].getWeight();
      if (egg == 0) {
        dp[egg][cap] = (size > cap) ? 0 : weight;
      } else {
        dp[egg][cap] = dp[egg - 1][cap];
        if (cap >= size)
          dp[egg][cap] =
              std::max(dp[egg][cap], dp[egg - 1][cap - size] + weight);
      }
    }
  }

  std::vector<Interval> getIntervals(int length, int pieces) {
    double interval = std::max(1., static_cast<double>(length) / pieces);
    int prev = -1, count = 0;
    double end = -1;
    std::vector<Interval> res;

    while (prev < length) {
      end += interval;
      count++;
      int next = std::min(static_cast<int>(end), length);
      if (count == pieces) next = length;

      res.push_back(std::make_pair(prev + 1, next));
      prev = next;
    }

    return res;
  }

  void mergeIntervals(std::vector<Interval>& intervals,
                      std::vector<GrainOfSand>& grains, int left, int right) {
    if (left == right) return;

    int mid = (left + right) / 2;
    std::future<void> future;
    future = councilOfShamans.enqueue([this, &intervals, &grains, left, mid]() {
      mergeIntervals(intervals, grains, left, mid);
    });
    mergeIntervals(intervals, grains, mid + 1, right);

    future.wait();
    merge(grains, intervals[left].first, intervals[mid + 1].first,
          intervals[right].second);
  }
};

#endif  // SRC_ADVENTURE_H_
