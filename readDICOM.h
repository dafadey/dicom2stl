#pragma once
#include <string>

#include <itkImage.h>
#include <itkGDCMSeriesFileNames.h>

#include "vgeo.h"

struct DICOMreader {
  using PixelType = signed short;
  using ImageType = itk::Image<PixelType, 3>;
  using NamesGeneratorType = itk::GDCMSeriesFileNames;

  DICOMreader(const std::string& _dirName, const std::string& _seriesId);

  bool read(std::vector<vgeo>& vgs);

  std::string dirName;
  std::string seriesId;
};
