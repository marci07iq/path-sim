#include "Main.h"

double qperm = 95788332; //1.6e-19 / 9.109e31 //95 788 332

double timestep = 1e-9;

double omega = qperm;

int symcount = 10;

int symlen = 1000;

int precision = 1;

dVec3 emitter = {0, 0, 0};

dVec3 emittererror = { 0, 0, 0 };

dVec3 velocity = {1e6, 0, 0};

float ktPerm = (300 * 1.381E-23) / (1.673E-27);

bool closeonend = false;

string Efilename = "Efield.txt";

string Bfilename = "Bfield.txt";

dVec3 Estep;
dVec3 Estart;
iVec3 Esize;
vector<vector<dVec3> > Efield;

dVec3 Bstep;
dVec3 Bstart;
iVec3 Bsize;
vector<vector<dVec3> > Bfield;

bool logpos = false;

ofstream logfile("log.txt");

int keepgoing = 0;

unsigned seed;
std::default_random_engine generator;

std::normal_distribution<double> distribution;

float pfmodf(float _n, float _q) {
  return _n - _q*floor(_n/_q);
}

dVec3 interpolate2D(vector<vector<dVec3> > &field, dVec3 &step, dVec3 &start, iVec3 &size, dVec3 pos, dVec3 flipX, dVec3 flipZ) {
  int ix = floor((pos.x - start.x) / step.x);
  int iz = floor((pos.z - start.z) / step.z);

  if (0 <= ix && ix + 1 < size.x) {
    if (0 <= iz && iz + 1 < size.z) {
      float qx = pfmodf(pos.x - start.x, step.x);
      float qz = pfmodf(pos.z - start.z, step.z);

      return
        (field[ix + 1][iz + 1] * (qx)* (qz)+
        field[ix + 1][iz] * (qx)* (1 - qz) +
        field[ix][iz + 1] * (1 - qx) * (qz)+
        field[ix][iz] * (1 - qx) * (1 - qz));
    }
  }

  if(ix < 0) {
    dVec3 pres = interpolate2D(field, step, start, size, { 2*start.x - pos.x, pos.y, pos.z }, flipX, flipZ);
    return {pres.x * flipX.x, pres.y* flipX.y, pres.z* flipX.z };
  }
  if (iz < 0) {
    dVec3 pres = interpolate2D(field, step, start, size, { pos.x, pos.y, 2*start.z-pos.z }, flipX, flipZ);
    return{ pres.x* flipZ.x, pres.y* flipZ.y, pres.z * flipZ.z };
  }

  if (ix + 1 >= size.x) {
    keepgoing = 1;
  }
  if (iz + 1 >= size.z) {
    keepgoing = 2;
  }

  //throw 1; //Out of bounds
  return {0,0,0};
}

dVec3 B(double time, dVec3 pos) {
  dVec3 norot =  interpolate2D(Bfield, Bstep, Bstart, Bsize, {sqrt(pos.x*pos.x + pos.y*pos.y), 0, pos.z}, {0, 0, 0}, {-1, -1, 1});
  float phi = atan2(pos.y, pos.x);
  return {cos(phi) * norot.x, sin(phi) * norot.x, norot.z};
  //return {0,0,1};
}

dVec3 E(double time, dVec3 pos) {
  int okg = keepgoing;
  dVec3 field = interpolate2D(Efield, Estep, Estart, Esize, {pos.x, 0, pos.z}, {1, 1, -1}, {1, 1, -1}) * 2e3 * cos(omega * time);
  if(keepgoing == 1) {keepgoing= okg;
  }
  return {field.x, field.y, field.z};
}

dVec3 acc(double time, dVec3 pos, dVec3 vel) {
  return (crs(vel, B(time, pos) )+  E(time, pos))*qperm; //(q v x B + q E)/mS
}

void rundekutta(double time, movement& path) {
  dVec3 dx1 = path.vel * timestep;
  dVec3 dv1 = acc(time, path.pos, path.vel) * timestep;
  dVec3 dx2 = (path.vel + dv1 / 2) * timestep;
  dVec3 dv2 = acc(time + timestep/2, path.pos + dx1 / 2, path.vel + dv1/2) * timestep;
  dVec3 dx3 = (path.vel + dv2 / 2) * timestep;
  dVec3 dv3 = acc(time + timestep / 2, path.pos + dx2 / 2, path.vel + dv2 / 2) * timestep;
  dVec3 dx4 = (path.vel + dv3) * timestep;
  dVec3 dv4 = acc(time + timestep, path.pos + dx3, path.vel + dv3) * timestep;
  dVec3 dx = (dx1 + dx2*2 + dx3*2 + dx4)/6;
  dVec3 dv = (dv1 + dv2 * 2 + dv3 * 2 + dv4) / 6;
  path = { path.pos + dx, path.vel + dv, acc(time, path.pos + dx, path.vel + dv), time + timestep};
}

void parsearg(string &s, stringstream &args) {
  if (s == "timestep") {
    args >> timestep;
    return;
  }
  if (s == "qperm") {
    args >> qperm;
    return;
  }
  if (s == "omega") {
    args >> omega;
    return;
  }
  if (s == "symcount") {
    args >> symcount;
    return;
  }
  if (s == "precision") {
    args >> precision;
    return;
  }
  if (s == "symlen") {
    args >> symlen;
    return;
  }
  if (s == "emitter") {
    args >> emitter;
    return;
  }
  if (s == "emittererror") {
    args >> emittererror;
    return;
  }
  if (s == "velocity") {
    args >> velocity;
    return;
  }
  if (s == "ktperm") {
    args >> ktPerm;
    return;
  }
  if (s == "closeonend") {
    args >> closeonend;
    return;
  }
  if (s == "Efilename") {
    args >> Efilename;
    return;
  }
  if (s == "Estep") {
    args >> Estep;
    return;
  }
  if (s == "Estart") {
    args >> Estart;
    return;
  }
  if (s == "Esize") {
    args >> Esize;
    return;
  }
  if (s == "Bfilename") {
    args >> Bfilename;
    return;
  }
  if (s == "Bstep") {
    args >> Bstep;
    return;
  }
  if (s == "Bstart") {
    args >> Bstart;
    return;
  }
  if (s == "Bsize") {
    args >> Bsize;
    return;
  }
  if (s == "logpos") {
    args >> logpos;
    return;
  }
  else {
    cout << "Unknown argument " << s << endl;
    logfile << "Unknown argument " << s << endl;
  }
}

double getPm1() {
  return 2*((generator() - generator.min()) / (generator.max()-generator.min()))-1;
}

double getGauss() {
  return distribution(generator);
}

void getparams() {
  ifstream in("in.txt");

  while (in.good()) {
    string s;
    getline(in, s);

    if (s.length() && s[0] != '#') { //not comment
      stringstream splitter;
      splitter << s;
      string argname;
      splitter >> argname;
      parsearg(argname, splitter);
    }
  }
}

void loadFields() {
  cout << "Loading E field from " << Efilename << endl;

  Efield = vector<vector<dVec3> > (Esize.x, vector<dVec3> (Esize.z));

  ifstream Ein(Efilename);

  while (Ein.good()) {
    string s;
    getline(Ein, s);

    if (s.length() && s[0] != '%') { //not comment
      stringstream splitter;
      splitter << s;
      double x, z, V, Ex, Ez;
      splitter >> x >> z >> V >> Ex >> Ez;
      Efield[round((x - Estart.x) / Estep.x)][round((z - Estart.z) / Estep.z)] = {Ex, 0, Ez};
    }
  }

  cout << "Loading B field from " << Bfilename << endl;

  Bfield = vector<vector<dVec3> >(Bsize.x, vector<dVec3>(Bsize.z));

  ifstream Bin(Bfilename);

  while (Bin.good()) {
    string s;
    getline(Bin, s);

    if (s.length() && s[0] != '%') { //not comment
      stringstream splitter;
      splitter << s;
      double x, z, Bx, Bz;
      splitter >> x >> z >> Bx >> Bz;
      Bfield[round((x - Bstart.x) / Bstep.x)][round((z - Bstart.z) / Bstep.z)] = { Bx, 0, Bz };
    }
  }
}

int main() {
  srand(time(0));

  seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator = default_random_engine(seed);

  distribution = normal_distribution<double>(0.0, sqrt(ktPerm));

  getparams();

  loadFields();

  movement path;

  ofstream of;
  if(logpos) {
    of = ofstream("out.txt");
  }
  ofstream collected("coll.txt");
  for(int n=0; n<symcount; n++) {
    path = {emitter + dVec3{getPm1() * emittererror.x, getPm1() * emittererror.y, getPm1() * emittererror.z}, velocity + dVec3{ abs(getGauss()) , getGauss() , getGauss() }, {0,0,0}, 0};
    
    cout << "Simulation " << n + 1 << " from " << symcount << endl;
    cout << "Time " << n / double(2*symcount * omega) << " Position " << path.pos << " Velocity " << path.vel << endl;
    logfile << "Simulation " << n + 1 << " from " << symcount << endl;
    logfile << "Time " << n / double(2 * symcount * omega) << " Position " << path.pos << " Velocity " << path.vel << endl;
    collected << n / double(2*symcount * omega) << "\t" << path.pos << "\t" << path.vel << "\t";


    keepgoing = 0;

    int i;

    for (i = 0; i < symlen && keepgoing == 0 && path.pos.x*path.pos.x + path.pos.y*path.pos.y < 0.01; i++) {
      if(logpos && i % precision == 0) {
        of << n << "\t" << path.pos << "\t" << path.t << endl;
      }

      rundekutta(i*timestep + n / double(2*symcount * omega), path);
      if (keepgoing == 1) {
        cout << "Particle exited accelerator." << endl;
        logfile << "Particle exited accelerator." << endl;
      }
      if (keepgoing == 2) {
        cout << "Particle collided with accelerator." << endl;
        logfile << "Particle collided with accelerator." << endl;
      }
    }

    collected << i*timestep + n / double(2*symcount * omega) << "\t" << path.pos << "\t" << path.vel << "\t" << keepgoing << endl;
  }

  cout << "Done" << endl;
  logfile << "Done" << endl;

  logfile.close();
  collected.close();
  if (logpos) {
    of.close();
  }

  if (!closeonend) {
    int n;
    cin >> n;
  }
  return 0;
}