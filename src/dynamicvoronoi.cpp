#include "dynamicvoronoi.h"

#include <math.h>
#include <iostream>

DynamicVoronoi::DynamicVoronoi() {
  sqrt2_ = sqrt(2.0);
  data_ = NULL;
  gridMap_ = NULL;
  alternativeDiagram_ = NULL;
  allocatedGridMap_ = false;
}

DynamicVoronoi::~DynamicVoronoi() {
  if (data_) {
    for (int x=0; x<sizeX_; x++) delete[] data_[x];
    delete[] data_;
  }
  if (allocatedGridMap_ && gridMap_) {
    for (int x=0; x<sizeX_; x++) delete[] gridMap_[x];
    delete[] gridMap_;
  }
}

void DynamicVoronoi::initializeEmpty(int _sizeX, int _sizeY, bool initGridMap) {
  //先清空历史数据
  if (data_) {
    for (int x=0; x<sizeX_; x++) delete[] data_[x];
    delete[] data_;
    data_ = NULL;
  }
  if(alternativeDiagram_){
    for (int x=0; x<sizeX_; x++) delete[] alternativeDiagram_[x];
    delete[] alternativeDiagram_;
    alternativeDiagram_ = NULL;
  }
  if (initGridMap) {
    if (allocatedGridMap_ && gridMap_) {
      for (int x=0; x<sizeX_; x++) delete[] gridMap_[x];
      delete[] gridMap_;
      gridMap_ = NULL;
      allocatedGridMap_ = false;
    }
  }

  //为数组分配内存
  sizeX_ = _sizeX;
  sizeY_ = _sizeY;
  data_ = new dataCell*[sizeX_];
  for (int x=0; x<sizeX_; x++) data_[x] = new dataCell[sizeY_];

  if (initGridMap) {
    gridMap_ = new bool*[sizeX_];
    for (int x=0; x<sizeX_; x++) gridMap_[x] = new bool[sizeY_];
    allocatedGridMap_ = true;
  }

  dataCell c;
  c.dist = INFINITY;
  c.sqdist = INT_MAX;
  c.obstX = invalidObstData;
  c.obstY = invalidObstData;
  c.voronoi = free;
  c.queueing = fwNotQueued;
  c.needsRaise = false;

  //为数组填充初始化数据
  for (int x=0; x<sizeX_; x++)
    for (int y=0; y<sizeY_; y++) data_[x][y] = c;

  if (initGridMap) {
    for (int x=0; x<sizeX_; x++)
      for (int y=0; y<sizeY_; y++) gridMap_[x][y] = 0;
  }
}

//输入二值地图gridmap，根据元素是否被占用，更新data_
void DynamicVoronoi::initializeMap(int _sizeX, int _sizeY, bool** _gridMap) {
  gridMap_ = _gridMap;
  initializeEmpty(_sizeX, _sizeY, false);

  for (int x=0; x<sizeX_; x++) {
    for (int y=0; y<sizeY_; y++) {
      if (gridMap_[x][y]) {             //如果gridmap_中的(x,y)被占用了
        dataCell c = data_[x][y];
        if (!isOccupied(x,y,c)) {       //如果c没有被占用，即data_中的(x,y)没被占用，需要更新
          bool isSurrounded = true;     //如果在gridmap_中的邻居元素全被占用，isSurrounded = true
          for (int dx=-1; dx<=1; dx++) {
            int nx = x+dx;
            if (nx<=0 || nx>=sizeX_-1) continue;
            for (int dy=-1; dy<=1; dy++) {
              if (dx==0 && dy==0) continue;
              int ny = y+dy;
              if (ny<=0 || ny>=sizeY_-1) continue;

              if (!gridMap_[nx][ny]) {  //如果在gridmap_中的邻居元素有任意一个没被占用（就是障碍物边界点），isSurrounded = false
                isSurrounded = false;
                break;
              }
            }
          }
          if (isSurrounded) {           //如果九宫格全部被占用
            c.obstX = x;
            c.obstY = y;
            c.sqdist = 0;
            c.dist=0;
            c.voronoi=occupied;
            c.queueing = fwProcessed;
            data_[x][y] = c;
          } else {
            setObstacle(x,y);       //不同之处在于：将(x,y)加入addList_
          }
        }
      }
    }
  }
}

//要同时更新gridmap和data_
void DynamicVoronoi::occupyCell(int x, int y) {
  gridMap_[x][y] = 1;     //更新gridmap
  setObstacle(x,y);
}

//要同时更新gridmap和data_
void DynamicVoronoi::clearCell(int x, int y) {
  gridMap_[x][y] = 0;     //更新gridmap
  removeObstacle(x,y);
}

//只更新data_
void DynamicVoronoi::setObstacle(int x, int y) {
  dataCell c = data_[x][y];
  if(isOccupied(x,y,c)) {                 //如果data_中的(x,y)被占用
    return;
  }

  addList_.push_back(INTPOINT(x,y));      //加入addList_
  c.obstX = x;
  c.obstY = y;
  data_[x][y] = c;
}

//只更新data_
void DynamicVoronoi::removeObstacle(int x, int y) {
  dataCell c = data_[x][y];
  if(isOccupied(x,y,c) == false) {          //如果data_中的(x,y)没有被占用，无需处理
    return;
  }

  removeList_.push_back(INTPOINT(x,y));     //将(x,y)加入removeList_
  c.obstX = invalidObstData;
  c.obstY  = invalidObstData;
  c.queueing = bwQueued;
  data_[x][y] = c;
}

//用新的障碍物信息替换旧的障碍物信息
//如果points为空，就是清除障碍物；
//初始时lastObstacles_为空，第一次调用exchangeObstacles()就是纯粹的添加障碍物
void DynamicVoronoi::exchangeObstacles(std::vector<INTPOINT>& points) {
  for (unsigned int i=0; i<lastObstacles_.size(); i++) {
    int x = lastObstacles_[i].x;
    int y = lastObstacles_[i].y;
    bool v = gridMap_[x][y];
    if (v) {    //如果(x,y)被占用了，不处理，怀疑这里逻辑反了。要移除旧的障碍物，这里应该是(!v)表示没被占用就不处理，占用了就移除
      continue;
    }
    removeObstacle(x,y);
  }
  lastObstacles_.clear();

  for (unsigned int i=0; i<points.size(); i++) {
    int x = points[i].x;
    int y = points[i].y;
    bool v = gridMap_[x][y];
    if (v) {    //如果(x,y)被占用了，不处理。否则，添加占用
      continue;
    }
    setObstacle(x,y);
    lastObstacles_.push_back(points[i]);
  }
}

void DynamicVoronoi::update(bool updateRealDist) {
  //将发生状态变化（占用<-->不占用）的元素加入open_优先队列
  commitAndColorize(updateRealDist);

  while (!open_.empty()) {
    INTPOINT p = open_.pop();
    int x = p.x;
    int y = p.y;
    dataCell c = data_[x][y];

    if(c.queueing==fwProcessed) {
      continue;
    }

    if (c.needsRaise) {
      // RAISE
      //2层for循环，考察8个邻居栅格
      for (int dx=-1; dx<=1; dx++) {
        int nx = x+dx;
        if (nx<=0 || nx>=sizeX_-1) continue;
        for (int dy=-1; dy<=1; dy++) {
          if (dx==0 && dy==0) continue;
          int ny = y+dy;
          if (ny<=0 || ny>=sizeY_-1) continue;
          dataCell nc = data_[nx][ny];
          //nc有最近障碍物 且 不raise
          if (nc.obstX!=invalidObstData && !nc.needsRaise) {
            //如果nc原来的最近障碍物消失了
            if(!isOccupied(nc.obstX, nc.obstY, data_[nc.obstX][nc.obstY])) {
              open_.push(nc.sqdist, INTPOINT(nx,ny));
              nc.queueing = fwQueued;     //fwQueued表示刚加入open_排队？
              nc.needsRaise = true;       //需要raise，并清理掉原来的最近障碍物信息
              nc.obstX = invalidObstData;
              nc.obstY = invalidObstData;
              if (updateRealDist) {
                nc.dist = INFINITY;
              }
              nc.sqdist = INT_MAX;
              data_[nx][ny] = nc;
            } else {        //如果nc原来的最近障碍物还存在
              if(nc.queueing != fwQueued){    //??
                open_.push(nc.sqdist, INTPOINT(nx,ny));
                nc.queueing = fwQueued;
                data_[nx][ny] = nc;
              }
            }
          }
        }
      }
      c.needsRaise = false;
      c.queueing = bwProcessed;     //bwProcessed表示8个邻居元素raise处理完毕？
      data_[x][y] = c;
    }
    else if (c.obstX != invalidObstData && isOccupied(c.obstX, c.obstY, data_[c.obstX][c.obstY])) {
      //c是被占据的
      // LOWER
      c.queueing = fwProcessed;     //fwProcessed表示8个邻居元素lower处理完毕？
      c.voronoi = occupied;

      for (int dx=-1; dx<=1; dx++) {
        int nx = x+dx;
        if (nx<=0 || nx>=sizeX_-1) continue;
        for (int dy=-1; dy<=1; dy++) {
          if (dx==0 && dy==0) continue;
          int ny = y+dy;
          if (ny<=0 || ny>=sizeY_-1) continue;
          dataCell nc = data_[nx][ny];
          if(!nc.needsRaise) {
            int distx = nx-c.obstX;
            int disty = ny-c.obstY;
            int newSqDistance = distx*distx + disty*disty;
            bool overwrite =  (newSqDistance < nc.sqdist);    //nc到c的最近障碍物 比 nc到其最近障碍物 更近
            if(!overwrite && newSqDistance==nc.sqdist) {
              //如果nc没有最近障碍物，或者 nc的最近障碍物消失了
              if (nc.obstX == invalidObstData || isOccupied(nc.obstX, nc.obstY, data_[nc.obstX][nc.obstY]) == false) {
                overwrite = true;
              }
            }
            if (overwrite) {
              open_.push(newSqDistance, INTPOINT(nx,ny));
              nc.queueing = fwQueued;     //fwQueued表示加入到open_等待lower()？
              if (updateRealDist) {
                nc.dist = sqrt((double) newSqDistance);
              }
              nc.sqdist = newSqDistance;
              nc.obstX = c.obstX;         //nc的最近障碍物 赋值为c的最近障碍物
              nc.obstY = c.obstY;
            } else {
              checkVoro(x,y,nx,ny,c,nc);
            }
            data_[nx][ny] = nc;
          }
        }
      }
    }
    data_[x][y] = c;
  }
}

float DynamicVoronoi::getDistance( int x, int y ) {
  if( (x>0) && (x<sizeX_) && (y>0) && (y<sizeY_)) {
    return data_[x][y].dist;
  }
  else return -INFINITY;
}

bool DynamicVoronoi::isVoronoi( int x, int y ) {
  dataCell c = data_[x][y];
  return (c.voronoi == free || c.voronoi == voronoiKeep);
}

bool DynamicVoronoi::isVoronoiAlternative(int x, int y) {
  int v = alternativeDiagram_[x][y];
  return (v == free || v == voronoiKeep);
}

void DynamicVoronoi::commitAndColorize(bool updateRealDist) {
  //addList_和removeList_中是触发Voronoi更新的元素，因此都要加入open_
  // ADD NEW OBSTACLES
  //addList_中都是障碍物边界点
  for (unsigned int i=0; i<addList_.size(); i++) {
    INTPOINT p = addList_[i];
    int x = p.x;
    int y = p.y;
    dataCell c = data_[x][y];

    if(c.queueing != fwQueued){
      if (updateRealDist) {
        c.dist = 0;
      }
      c.sqdist = 0;
      c.obstX = x;
      c.obstY = y;
      c.queueing = fwQueued;          //已加入open_优先队列
      c.voronoi = occupied;
      data_[x][y] = c;
      open_.push(0, INTPOINT(x,y));   //加入open_优先队列，加入open_的都是要更新的
    }
  }

  // REMOVE OLD OBSTACLES
  //removeList_中是要清除的障碍物栅格
  for (unsigned int i=0; i<removeList_.size(); i++) {
    INTPOINT p = removeList_[i];
    int x = p.x;
    int y = p.y;
    dataCell c = data_[x][y];

    //removeList_中对应的元素在data_中已经更新过，解除了占用
    //如果这里又出现了该元素被占用，说明是后来加入的，这里不处理
    if (isOccupied(x,y,c) == true) {
      continue; // obstacle was removed and reinserted
    }
    open_.push(0, INTPOINT(x,y));   //加入open_优先队列
    if (updateRealDist) {
      c.dist  = INFINITY;
    }
    c.sqdist = INT_MAX;
    c.needsRaise = true;            //因为清除了障碍物，最近障碍物距离要更新-增加
    data_[x][y] = c;
  }
  removeList_.clear();
  addList_.clear();
}


void DynamicVoronoi::checkVoro(int x, int y, int nx, int ny, dataCell& c, dataCell& nc) {
  if ((c.sqdist>1 || nc.sqdist>1) && nc.obstX!=invalidObstData) {
    if (abs(c.obstX-nc.obstX) > 1 || abs(c.obstY-nc.obstY) > 1) {
      //compute dist from x,y to obstacle of nx,ny
      int dxy_x = x-nc.obstX;
      int dxy_y = y-nc.obstY;
      int sqdxy = dxy_x*dxy_x + dxy_y*dxy_y;
      int stability_xy = sqdxy - c.sqdist;
      if (sqdxy - c.sqdist<0) return;

      //compute dist from nx,ny to obstacle of x,y
      int dnxy_x = nx - c.obstX;
      int dnxy_y = ny - c.obstY;
      int sqdnxy = dnxy_x*dnxy_x + dnxy_y*dnxy_y;
      int stability_nxy = sqdnxy - nc.sqdist;
      if (sqdnxy - nc.sqdist <0) return;

      //which cell is added to the Voronoi diagram?
      if(stability_xy <= stability_nxy && c.sqdist>2) {
        if (c.voronoi != free) {
          c.voronoi = free;
          reviveVoroNeighbors(x,y);
          pruneQueue_.push(INTPOINT(x,y));
        }
      }
      if(stability_nxy <= stability_xy && nc.sqdist>2) {
        if (nc.voronoi != free) {
          nc.voronoi = free;
          reviveVoroNeighbors(nx,ny);
          pruneQueue_.push(INTPOINT(nx,ny));
        }
      }
    }
  }
}


void DynamicVoronoi::reviveVoroNeighbors(int &x, int &y) {
  for (int dx=-1; dx<=1; dx++) {
    int nx = x+dx;
    if (nx<=0 || nx>=sizeX_-1) continue;
    for (int dy=-1; dy<=1; dy++) {
      if (dx==0 && dy==0) continue;
      int ny = y+dy;
      if (ny<=0 || ny>=sizeY_-1) continue;
      dataCell nc = data_[nx][ny];
      if (nc.sqdist != INT_MAX && !nc.needsRaise && (nc.voronoi == voronoiKeep || nc.voronoi == voronoiPrune)) {
        nc.voronoi = free;
        data_[nx][ny] = nc;
        pruneQueue_.push(INTPOINT(nx,ny));
      }
    }
  }
}


bool DynamicVoronoi::isOccupied(int x, int y) {
  dataCell c = data_[x][y];
  return (c.obstX==x && c.obstY==y);
}

bool DynamicVoronoi::isOccupied(int &x, int &y, dataCell &c) {
  return (c.obstX==x && c.obstY==y);
}

void DynamicVoronoi::visualize(const char *filename) {
  // write ppm files

  FILE* F = fopen(filename, "w");
  if (!F) {
    std::cerr << "could not open 'result.ppm' for writing!\n";
    return;
  }
  fprintf(F, "P6\n#\n");
  fprintf(F, "%d %d\n255\n", sizeX_, sizeY_);

  //fputc()执行3次，其实是依次对一个像素的RGB颜色赋值
  for(int y = sizeY_-1; y >=0; y--){
    for(int x = 0; x<sizeX_; x++){
      unsigned char c = 0;
      if (alternativeDiagram_!=NULL && (alternativeDiagram_[x][y] == free || alternativeDiagram_[x][y]==voronoiKeep)) {
        //和alternative模式相关，先不用管
        fputc( 255, F );
        fputc( 0, F );
        fputc( 0, F );
      } else if(isVoronoi(x,y)){  //画Voronoi边
        fputc( 0, F );
        fputc( 0, F );
        fputc( 255, F );
      } else if (data_[x][y].sqdist==0) {  //填充障碍物
        fputc( 0, F );
        fputc( 0, F );
        fputc( 0, F );
      } else {    //填充Voronoi区块内部
        float f = 80+(sqrt(data_[x][y].sqdist)*10);
        if (f>255) f=255;
        if (f<0) f=0;
        c = (unsigned char)f;
        fputc( c, F );
        fputc( c, F );
        fputc( c, F );
      }
    }
  }
  fclose(F);
}


void DynamicVoronoi::prune() {
  // filler
  //先遍历pruneQueue_中的元素，判断是否要加入到sortedPruneQueue_，（为什么要这一步？？？）
  //再遍历sortedPruneQueue_中的元素，判断其是剪枝、保留、重试。
  while(!pruneQueue_.empty()) {
    INTPOINT p = pruneQueue_.front();
    pruneQueue_.pop();
    int x = p.x;
    int y = p.y;

    if (data_[x][y].voronoi==occupied) continue;    //如果(x,y)是occupied，无需处理，不可能是Voronoi
    if (data_[x][y].voronoi==freeQueued) continue;  //如果(x,y)是freeQueued，已经加入到sortedPruneQueue_，略过

    data_[x][y].voronoi = freeQueued;
    sortedPruneQueue_.push(data_[x][y].sqdist, p);

    /* tl t tr
       l c r
       bl b br */

    dataCell tr,tl,br,bl;
    tr = data_[x+1][y+1];
    tl = data_[x-1][y+1];
    br = data_[x+1][y-1];
    bl = data_[x-1][y-1];

    dataCell r,b,t,l;
    r = data_[x+1][y];
    l = data_[x-1][y];
    t = data_[x][y+1];
    b = data_[x][y-1];

    if (x+2<sizeX_ && r.voronoi==occupied) {
      // fill to the right
      //如果r的上下左右4个元素都!=occupied
      if (tr.voronoi!=occupied && br.voronoi!=occupied && data_[x+2][y].voronoi!=occupied) {
        r.voronoi = freeQueued;
        sortedPruneQueue_.push(r.sqdist, INTPOINT(x+1,y));
        data_[x+1][y] = r;
      }
    }
    if (x-2>=0 && l.voronoi==occupied) {
      // fill to the left
      //如果l的上下左右4个元素都!=occupied
      if (tl.voronoi!=occupied && bl.voronoi!=occupied && data_[x-2][y].voronoi!=occupied) {
        l.voronoi = freeQueued;
        sortedPruneQueue_.push(l.sqdist, INTPOINT(x-1,y));
        data_[x-1][y] = l;
      }
    }
    if (y+2<sizeY_ && t.voronoi==occupied) {
      // fill to the top
      //如果t的上下左右4个元素都!=occupied
      if (tr.voronoi!=occupied && tl.voronoi!=occupied && data_[x][y+2].voronoi!=occupied) {
        t.voronoi = freeQueued;
        sortedPruneQueue_.push(t.sqdist, INTPOINT(x,y+1));
        data_[x][y+1] = t;
      }
    }
    if (y-2>=0 && b.voronoi==occupied) {
      // fill to the bottom
      //如果b的上下左右4个元素都!=occupied
      if (br.voronoi!=occupied && bl.voronoi!=occupied && data_[x][y-2].voronoi!=occupied) {
        b.voronoi = freeQueued;
        sortedPruneQueue_.push(b.sqdist, INTPOINT(x,y-1));
        data_[x][y-1] = b;
      }
    }
  }

  while(!sortedPruneQueue_.empty()) {
    INTPOINT p = sortedPruneQueue_.pop();
    dataCell c = data_[p.x][p.y];
    int v = c.voronoi;
    if (v!=freeQueued && v!=voronoiRetry) { // || v>free || v==voronoiPrune || v==voronoiKeep) {
      //      assert(v!=retry);
      continue;
    }

    markerMatchResult r = markerMatch(p.x,p.y);
    if (r==pruned) {
      c.voronoi = voronoiPrune;     //对(x,y)即c剪枝
    }
    else if (r==keep) {
      c.voronoi = voronoiKeep;      //对(x,y)即c保留，成为Voronoi的边
    }
    else {
      c.voronoi = voronoiRetry;
      pruneQueue_.push(p);
    }
    data_[p.x][p.y] = c;

    //把需要retry的元素由pruneQueue_转移到sortedPruneQueue_
    //这样可以继续本while()循环，直到pruneQueue_和sortedPruneQueue_都为空
    if (sortedPruneQueue_.empty()) {
      while (!pruneQueue_.empty()) {
        INTPOINT p = pruneQueue_.front();
        pruneQueue_.pop();
        sortedPruneQueue_.push(data_[p.x][p.y].sqdist, p);
      }
    }
  }
}

void DynamicVoronoi::updateAlternativePrunedDiagram() {

  if(alternativeDiagram_==NULL){
    alternativeDiagram_ = new int*[sizeX_];
    for(int x=0; x<sizeX_; x++){
      alternativeDiagram_[x] = new int[sizeY_];
    }
  }


  std::queue<INTPOINT> end_cells;
  BucketPrioQueue<INTPOINT> sortedPruneQueue;
  for(int x=1; x<sizeX_-1; x++){
    for(int y=1; y<sizeY_-1; y++){
      dataCell& c = data_[x][y];
	alternativeDiagram_[x][y] = c.voronoi;
	if(c.voronoi <=free){
	  sortedPruneQueue.push(c.sqdist, INTPOINT(x,y));
	  end_cells.push(INTPOINT(x, y));
	}
    }
  }

  for(int x=1; x<sizeX_-1; x++){
    for(int y=1; y<sizeY_-1; y++){
      if( getNumVoronoiNeighborsAlternative(x, y) >= 3){
	alternativeDiagram_[x][y] = voronoiKeep;
	sortedPruneQueue.push(data_[x][y].sqdist, INTPOINT(x,y));
	end_cells.push(INTPOINT(x, y));
      }
    }
  }

  for(int x=1; x<sizeX_-1; x++){
    for(int y=1; y<sizeY_-1; y++){
      if( getNumVoronoiNeighborsAlternative(x, y) >= 3){
	alternativeDiagram_[x][y] = voronoiKeep;
	sortedPruneQueue.push(data_[x][y].sqdist, INTPOINT(x,y));
	end_cells.push(INTPOINT(x, y));
      }
    }
  }


  while (!sortedPruneQueue.empty()) {
    INTPOINT p = sortedPruneQueue.pop();

    if (markerMatchAlternative(p.x, p.y)) {
      alternativeDiagram_[p.x][p.y]=voronoiPrune;
    } else {
  alternativeDiagram_[p.x][p.y]=voronoiKeep;
    }
  }

  // //delete worms
  while (!end_cells.empty()) {
    INTPOINT p = end_cells.front();
    end_cells.pop();

    if (isVoronoiAlternative(p.x,p.y) && getNumVoronoiNeighborsAlternative(p.x, p.y) == 1) {
      alternativeDiagram_[p.x][p.y] = voronoiPrune;

      for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
          if (!(dx || dy) || (dx && dy)) {
            continue;
          }
          int nx = p.x + dx;
          int ny = p.y + dy;
          if (nx < 0 || nx >= sizeX_ || ny < 0 || ny >= sizeY_) {
            continue;
          }
          if (isVoronoiAlternative(nx,ny)) {
            if (getNumVoronoiNeighborsAlternative(nx, ny) == 1) {
              end_cells.push(INTPOINT(nx, ny));
            }
          }
        }
      }
    }
  }
}

bool DynamicVoronoi::markerMatchAlternative(int x, int y) {
// prune if this returns true

  bool f[8];

  int nx, ny;
  int dx, dy;

  int i = 0;
//  int obstacleCount=0;
  int voroCount = 0;
  for (dy = 1; dy >= -1; dy--) {
    ny = y + dy;
    for (dx = -1; dx <= 1; dx++) {
      if (dx || dy) {
        nx = x + dx;
        int v = alternativeDiagram_[nx][ny];
        bool b = (v <= free && v != voronoiPrune);
        //	if (v==occupied) obstacleCount++;
        f[i] = b;
        if (v <= free && !(dx && dy))
          voroCount++;
        i++;
      }
    }
  }

  /*
   * 5 6 7
   * 3   4
   * 0 1 2
   */

  {
    //connected horizontal or vertically to only one cell
    if (voroCount == 1 && (f[1] || f[3] || f[4] || f[6])) {
      return false;
    }

    // 4-connected
    if ((!f[0] && f[1] && f[3]) || (!f[2] && f[1] && f[4]) || (!f[5] && f[3] && f[6]) || (!f[7] && f[6] && f[4]))
      return false;

    if ((f[3] && f[4] && !f[1] && !f[6]) || (f[1] && f[6] && !f[3] && !f[4]))
      return false;

  }
  return true;
}

int DynamicVoronoi::getNumVoronoiNeighborsAlternative(int x, int y) {
  int count = 0;
  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      if ((dx == 0 && dy == 0) || (dx != 0 && dy != 0)) {
        continue;
      }

      int nx = x + dx;
      int ny = y + dy;
      if (nx < 0 || nx >= sizeX_ || ny < 0 || ny >= sizeY_) {
        continue;
      }
      if (alternativeDiagram_[nx][ny]==free || alternativeDiagram_[nx][ny]==voronoiKeep) {
        count++;
      }
    }
  }
  return count;
}

//根据(x,y)邻居栅格的连接模式，判断是否要对(x,y)剪枝
DynamicVoronoi::markerMatchResult DynamicVoronoi::markerMatch(int x, int y) {
  // implementation of connectivity patterns
  bool f[8];
  int nx, ny;
  int dx, dy;
  int i=0;
  //voroCount是对所有邻居栅格的统计，voroCountFour是对上下左右4个邻居栅格的统计
  int voroCount=0;
  int voroCountFour=0;

  for (dy=1; dy>=-1; dy--) {
    ny = y+dy;
    for (dx=-1; dx<=1; dx++) {
      if (dx || dy) {   //不考虑(x,y)点
        nx = x+dx;
        dataCell nc = data_[nx][ny];
        int v = nc.voronoi;
        bool b = (v<=free && v!=voronoiPrune);    //既不是occupied又不是voronoiPrune，即可能保留的栅格
        f[i] = b;
        if (b) {
          voroCount++;
          if (!(dx && dy)) {      //对上下左右4个点
            voroCountFour++;
          }
        }
        i++;
      }
    }
  }
  // i和位置的对应关系如下：
  //    | 0 | 1 | 2 |
  //    | 3 |   | 4 |
  //    | 5 | 6 | 7 |
  //8个邻居栅格中最多有2个，上下左右只有1个可能保留的栅格
  if (voroCount<3 && voroCountFour==1 && (f[1] || f[3] || f[4] || f[6])) {
    return keep;
  }

  // 4-connected
  //    | 0 | 1 | ? |               | ? | 1 | 0 |             | ? | ? | ? |             | ? | ? | ? |
  //    | 1 |   | ? |               | ? |   | 1 |             | 1 |   | ? |             | ? |   | 1 |
  //    | ? | ? | ? |               | ? | ? | ? |             | 0 | 1 | ? |             | ? | 1 | 0 |
  //对应《Efficient Grid-Based Spatial Representations for Robot Navigation in Dynamic Environments》中的4-connected P14模式，旋转3次90度
  if ((!f[0] && f[1] && f[3]) || (!f[2] && f[1] && f[4]) || (!f[5] && f[3] && f[6]) || (!f[7] && f[6] && f[4])) return keep;

  //    | ? | 0 | ? |                       | ? | 1 | ? |
  //    | 1 |   | 1 |                       | 0 |   | 0 |
  //    | ? | 0 | ? |                       | ? | 1 | ? |
  //对应文章中的4-connected P24模式，旋转1次90度
  if ((f[3] && f[4] && !f[1] && !f[6]) || (f[1] && f[6] && !f[3] && !f[4])) return keep;

  // keep voro cells inside of blocks and retry later
  //(x,y)周围可能保留的栅格很多，此时无法判断是否要对(x,y)剪枝
  if (voroCount>=5 && voroCountFour>=3 && data_[x][y].voronoi!=voronoiRetry) {
    return retry;
  }

  return pruned;
}
