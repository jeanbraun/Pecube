"""
    This file gather all functions and classes related to the GIS.

"""

from osgeo import gdal
from osgeo import osr, ogr
import rioxarray
import xarray as xr
import numpy as np
import math

#----------------------------------------------------------
def load_raster_file(path):
    """
    This function enables to load raster files.
    It returns elevation, and coordinates.

    """
    # # Temporary copy the raster file
    
    # with rasterio.open(str(path)) as dataFile:
    #     data = dataFile.load()
    #     print(data)
    # # Remove temp file and clean
    # os.remove(tempFile)
    # # Check whether the raster file is in geographic coordinate (WGS84)
    # try:
    #     data_reproject = data.rio.reproject("EPSG:4326")
    #     print(data_reproject)
    # except:
    #     QErrorMessage().showMessage("The DEM file cannot be read.")
    #     return
    data = gdal.Open(path,gdal.GA_ReadOnly)
    ProjSys = data.GetProjection()
    srs = osr.SpatialReference(wkt=ProjSys)
    CoordSys = srs.GetAttrValue('AUTHORITY',1)
    
    # Get size of the DEM
    nx = data.RasterXSize
    ny = data.RasterYSize
    geotransform = data.GetGeoTransform()
    lon0, dlon, dxdy, lat0, dydx, dlat = data.GetGeoTransform()
    projection = data.GetProjection()
    
    #coordinates of the bottom left corner
    lat0 = lat0 +(data.RasterYSize * dlat)
    if srs.GetAttrValue('AUTHORITY',1) != '4326':
        xtemp = [lon0, lon0+dlon]; ytemp = [lat0, lat0+dlat]
        newCoord = TransformCoord(int(srs.GetAttrValue('AUTHORITY',1)),4326,xtemp,ytemp)
        dlat = newCoord.newCoord[0][1] - newCoord.newCoord[1][1]
        dlon = newCoord.newCoord[0][0] - newCoord.newCoord[1][0]
        lon0 = newCoord.newCoord[0][0]; lat0 = newCoord.newCoord[0][1]
    dlon = abs(dlon)
    dlat = abs(dlat)

    # Get elevation
    band = data.GetRasterBand(1)
    xsize = band.XSize
    ysize = band.YSize
    elevation_temp = band.ReadAsArray()
    elevation = np.reshape(np.flipud(elevation_temp),(nx*ny,1))
    # Handle nan values
    elevation[np.isnan(elevation)] = np.min(elevation[:])
    elevation[elevation<-5000] = 0
    data = None
    band = None
    
    return elevation,nx,ny,dlon,dlat,lon0,lat0


#----------------------------------------------------------
def km2lat(km):
    """to convert a distance in km to latitude in degree."""
    # Rayon de la Terre en km
    R = 6371.0
    # 1 degré de latitude correspond à environ 111 km
    km_per_degree_latitude = 111.11
    # Conversion
    degrees = km / km_per_degree_latitude
    return degrees

#----------------------------------------------------------
def km2lon(km, latitude):
    """to convert a distance in km to longitude in degree."""
    # Rayon de la Terre en km
    R = 6371.0
    # Conversion de la latitude en radians
    lat_rad = math.radians(latitude)
    # 1 degré de longitude en km à une latitude donnée
    km_per_degree_longitude = 111.11 * math.cos(lat_rad)
    # Conversion
    degrees = km / km_per_degree_longitude
    return degrees
    


#############################################################
class TransformCoord(object):
    """ This class take a Raster file with any projection system and convert it to geographic system WGS84"""
    def __init__(self,source,target,xcoord,ycoord):
        self.source = osr.SpatialReference()
        self.target = osr.SpatialReference()
        self.source.ImportFromEPSG(source)
        self.target.ImportFromEPSG(target)
        #self.coordTransform = osr.CoordinateTransformation(self.source,self.target)
        self.newCoord = []
        for i in range(len(xcoord)):
            self.newCoord.append(self.transform(xcoord[i],ycoord[i]))
        
    def transform(self,xcoord,ycoord):
        
        geom = ogr.Geometry(ogr.wkbPoint)
        geom.SetPoint_2D(0,xcoord,ycoord)
        geom.AssignSpatialReference(self.source)
        geom.TransformTo(self.target)
        return geom.GetPoint_2D()