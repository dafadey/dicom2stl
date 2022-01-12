#include "readDICOM.h"

#include <itkGDCMImageIO.h>
#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itkNrrdImageIOFactory.h>
#include <itkPNGImageIOFactory.h>
#include <itkVTKImageIOFactory.h>

DICOMreader::DICOMreader(const std::string& _dirName, const std::string& _seriesId) : dirName(_dirName), seriesId(_seriesId) {
  itk::NrrdImageIOFactory::RegisterOneFactory();
  itk::PNGImageIOFactory::RegisterOneFactory();
  itk::VTKImageIOFactory::RegisterOneFactory();
}

bool DICOMreader::read(std::vector<vgeo>& vgs) {
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails(true);
  nameGenerator->AddSeriesRestriction("0008|0021");
  nameGenerator->SetGlobalWarningDisplay(false);
  nameGenerator->SetDirectory(dirName.c_str());

  const std::vector<std::string>& seriesUID = nameGenerator->GetSeriesUIDs();
  auto seriesItr = seriesUID.begin();
  auto seriesEnd = seriesUID.end();

  if (seriesItr != seriesEnd)
  {
    std::cout << "The directory: ";
    std::cout << dirName << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl;
  }
  else
  {
    std::cout << "No DICOMs in: " << dirName << std::endl;
    return false;
  }

  while (seriesItr != seriesEnd)
  {
    std::cout << seriesItr->c_str() << std::endl;
    ++seriesItr;
  }

  seriesItr = seriesUID.begin();
  while (seriesItr != seriesUID.end())
  {
    std::string seriesIdentifier;
    if (!seriesId.empty()) // If seriesIdentifier given convert only that
    {
      seriesIdentifier = seriesId.c_str();
      seriesItr = seriesUID.end();
    }
    else // otherwise convert everything
    {
      seriesIdentifier = seriesItr->c_str();
      seriesItr++;
    }
    std::cout << "\nReading: ";
    std::cout << seriesIdentifier << std::endl;
    std::vector<std::string> fileNames = nameGenerator->GetFileNames(seriesIdentifier);

    using ReaderType = itk::ImageSeriesReader<ImageType>;
    ReaderType::Pointer reader = ReaderType::New();
    using ImageIOType = itk::GDCMImageIO;
    ImageIOType::Pointer dicomIO = ImageIOType::New();
    reader->SetImageIO(dicomIO);
    reader->SetFileNames(fileNames);
    reader->ForceOrthogonalDirectionOff(); // properly read CTs with gantry tilt

    ImageType* img = reader->GetOutput();
    img->Update();

    int nx = img->GetBufferedRegion().GetSize()[0];
    int ny = img->GetBufferedRegion().GetSize()[1];
    int nz = img->GetBufferedRegion().GetSize()[2];
    std::cout << "dim: " << nx << " x " << ny << " x " << nz << '\n';

    double dx = img->GetSpacing()[0];
    double dy = img->GetSpacing()[1];
    double dz = img->GetSpacing()[2];
    std::cout << "dd: " << dx << " x " << dy << " x " << dz << '\n';
    vgs.push_back(vgeo(nx,ny,nz));
    auto& vg = vgs[vgs.size()-1];
    vg.d[0] = dx;
    vg.d[1] = dy;
    vg.d[2] = dz;
    for(int i=0;i<3;i++)
      vg.o[i] = .5f * vg.d[i] * static_cast<float>(vg.dim[i]);
    for (int k = 0; k < nz; k++) {
      for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
          float v = static_cast<float>(img->GetPixel({ i,j,k }));
          vg[i+j*vg.dim[0]+k*vg.dim10] = v;
          vg.min = vg.min < v ? vg.min : v;
          vg.max = vg.max > v ? vg.max : v;
        }
      }
    }
  }
  return true;
}
