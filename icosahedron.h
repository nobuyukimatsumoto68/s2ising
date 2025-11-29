#pragma once


constexpr int NPatches = 10;
constexpr int PatchIdxS = 10;
constexpr int PatchIdxN = 11;
constexpr int PatchIdxVoid = 12;
constexpr int NIcosVertices = 12;

using IcosVertices = std::array<V3, NIcosVertices>;
using PatchEdges = std::array<std::array<int,4>, NPatches>;




struct Icosahedron {
  const double varphi = 0.5*( 1.0 + std::sqrt(5.0) );
  const double c1 = 0.5/varphi;
  const double s1 = 0.5*std::sqrt(std::sqrt(5.)*varphi);
  const double c2 = -0.5*varphi;
  const double s2 = 0.5*std::sqrt(std::sqrt(5.)/varphi);
  const double c3 = 1.0/std::sqrt(5.);
  const double s3 = 2.0/std::sqrt(5.);

  const double len = std::sqrt( 2.0-2.0/std::sqrt(5.0) );

  IcosVertices icos;
  PatchEdges patches;

  Icosahedron(){
    set_IcosVertices();
    set_PatchEdges();
  }

  // -----------------------------

  void set_IcosVertices(){
    icos[0] << 0.,0.,-1.;

    icos[1] << -s3, 0., -c3;
    icos[2] << -s3*c1, -s3*s1, -c3;
    icos[3] << -s3*c2, -s3*s2, -c3;
    icos[4] << -s3*c2, s3*s2, -c3;
    icos[5] << -s3*c1, s3*s1, -c3;

    icos[6] << s3*c2, -s3*s2, c3;
    icos[7] << s3*c1, -s3*s1, c3;
    icos[8] << s3, 0., c3;
    icos[9] << s3*c1, s3*s1, c3;
    icos[10] << s3*c2, s3*s2, c3;

    icos[11] << 0.,0.,1.;
  }

  void set_PatchEdges(){
    patches[0] = std::array<int,4>{1, 0, 6, 2};
    patches[1] = std::array<int,4>{2, 0, 7, 3};
    patches[2] = std::array<int,4>{3, 0, 8, 4};
    patches[3] = std::array<int,4>{4, 0, 9, 5};
    patches[4] = std::array<int,4>{5, 0, 10,1};

    patches[5] = std::array<int,4>{6, 2, 11, 7};
    patches[6] = std::array<int,4>{7, 3, 11, 8};
    patches[7] = std::array<int,4>{8, 4, 11, 9};
    patches[8] = std::array<int,4>{9, 5, 11, 10};
    patches[9] = std::array<int,4>{10,1, 11, 6};
  }

  V3 operator[](const Idx i) const { return icos[i]; }


};


