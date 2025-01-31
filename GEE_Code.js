// Advanced Evapotranspiration Analysis for Sudan
// Configuration and Constants
var CONFIG = {
  START_YEAR: 2000,
  END_YEAR: 2023,
  SCALE: 500,
  DROUGHT_THRESHOLDS: {
    EXTREME: -2,
    SEVERE: -1.5,
    MODERATE: -1
  }
};

// Color palettes
var PALETTES = {
  ET: [
    '#ffffcc', '#c7e9b4', '#7fcdbb', 
    '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84'
  ],
  DROUGHT: ['#d73027', '#fc8d59', '#fee090', '#e0f3f8']
};

// Load base datasets
var sudan = ee.FeatureCollection('FAO/GAUL/2015/level0')
    .filter(ee.Filter.eq('ADM0_NAME', 'Sudan'));
var admin1 = ee.FeatureCollection('FAO/GAUL/2015/level1')
    .filter(ee.Filter.eq('ADM0_NAME', 'Sudan'));

// Load and prepare MODIS ET data
var collection = ee.ImageCollection('MODIS/061/MOD16A2GF')
    .filterBounds(sudan);

// Calculate date range
var range = collection.reduceColumns(ee.Reducer.minMax(), ['system:time_start']);
var startDate = ee.Date(range.get('min'));
var endDate = ee.Date(range.get('max'));
var monthDiff = endDate.difference(startDate, 'months');

// Calculate monthly means with proper scaling
var calculateMonthlyET = function() {
  return ee.List.sequence(0, monthDiff).map(function(n) {
    var start = startDate.advance(n, 'month');
    var end = start.advance(1, 'month');
    
    return collection
      .filterDate(start, end)
      .mean()
      .multiply(0.1/8)  // Convert to mm/day
      .clip(sudan)
      .set({
        'system:time_start': start,
        'month': start.get('month'),
        'year': start.get('year')
      });
  });
};

var monthlyET = ee.ImageCollection(calculateMonthlyET());

// Statistical Analysis Functions
var calculateMannKendall = function(collection) {
  var sorted = collection.sort('system:time_start');
  var n = sorted.size();
  
  var kendall = sorted.map(function(image) {
    var time = ee.Number(image.get('system:time_start'));
    var laterImages = sorted.filterMetadata('system:time_start', 'greater_than', time);
    
    var increases = laterImages.map(function(laterImage) {
      return ee.Image(laterImage).subtract(image).gt(0);
    });
    
    return increases.sum();
  }).sum();
  
  var s = kendall.subtract(n.multiply(n.subtract(1)).divide(4));
  var variance = n.multiply(n.subtract(1)).multiply(n.multiply(2).add(5)).divide(72);
  var z = s.divide(variance.sqrt());
  
  return z;
};

// Drought Analysis
var calculateSETI = function(collection) {
  var mean = collection.mean();
  var stdDev = collection.reduce(ee.Reducer.stdDev());
  
  return collection.map(function(image) {
    var seti = image.subtract(mean)
                    .divide(stdDev)
                    .rename('SETI');
    
    var droughtClass = seti.expression(
      '(SETI <= extreme) ? 1 : (SETI <= severe) ? 2 : (SETI <= moderate) ? 3 : 4',
      {
        'SETI': seti.select('SETI'),
        'extreme': CONFIG.DROUGHT_THRESHOLDS.EXTREME,
        'severe': CONFIG.DROUGHT_THRESHOLDS.SEVERE,
        'moderate': CONFIG.DROUGHT_THRESHOLDS.MODERATE
      }
    ).rename('drought_class');
    
    return image.addBands([seti, droughtClass]);
  });
};

// Climate Data Integration
var addClimateData = function() {
  var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
      .filterBounds(sudan)
      .filterDate(startDate, endDate);
      
  var monthlyPrecip = ee.List.sequence(0, monthDiff).map(function(n) {
    var start = startDate.advance(n, 'month');
    var end = start.advance(1, 'month');
    
    return chirps
      .filterDate(start, end)
      .sum()
      .clip(sudan)
      .set('system:time_start', start);
  });
  
  return ee.ImageCollection(monthlyPrecip);
};

// Calculate all metrics
var ET = monthlyET.select('ET');
var SETI_collection = calculateSETI(ET);
var precipitation = addClimateData();
var trend = calculateMannKendall(ET);

// Visualization parameters
var etVis = {
  min: 0,
  max: 200,
  palette: PALETTES.ET
};

var droughtVis = {
  min: 1,
  max: 4,
  palette: PALETTES.DROUGHT
};

// Create UI elements
var createLegend = function() {
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '8px 15px'
    }
  });
  
  var title = ui.Label({
    value: 'ET (mm/day)',
    style: {
      fontWeight: 'bold',
      fontSize: '16px',
      margin: '0 0 4px 0'
    }
  });
  
  legend.add(title);
  
  var makeRow = function(color, name) {
    var colorBox = ui.Label({
      style: {
        backgroundColor: color,
        padding: '8px',
        margin: '0 0 4px 0'
      }
    });
    
    var description = ui.Label({
      value: name,
      style: {margin: '0 0 4px 6px'}
    });
    
    return ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.Flow('horizontal')
    });
  };
  
  var step = (etVis.max - etVis.min) / (PALETTES.ET.length - 1);
  PALETTES.ET.forEach(function(color, i) {
    var value = etVis.min + (step * i);
    legend.add(makeRow(color, value.toFixed(1)));
  });
  
  return legend;
};

// Create time series chart
var createTimeSeries = function() {
  return ui.Chart.image.series({
    imageCollection: ET,
    region: sudan,
    reducer: ee.Reducer.mean(),
    scale: CONFIG.SCALE,
    xProperty: 'system:time_start'
  })
  .setChartType('ScatterChart')
  .setOptions({
    title: 'Evapotranspiration in Sudan (MODIS 500m)',
    vAxis: {
      title: 'Average ET (mm/day)',
      viewWindow: {min: 0}
    },
    hAxis: {title: 'Date'},
    lineWidth: 1.5,
    pointSize: 2,
    trendlines: {0: {color: 'red'}},
    colors: ['#2c7fb8']
  });
};

// Calculate zonal statistics
var calculateZonalStats = function() {
  return ET.mean().reduceRegions({
    collection: admin1,
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    scale: CONFIG.SCALE
  });
};

// Export data
var exportData = function() {
  // 1. Export Administrative Level-1 Statistics
  Export.table.toDrive({
    collection: calculateZonalStats(),
    description: 'Sudan_ET_Admin1_Stats',
    fileFormat: 'CSV',
    folder: 'Sudan_ET_Analysis'
  });
  
  // 2. Export Monthly Time Series Data
  var timeSeriesData = ET.map(function(image) {
    var mean = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: sudan,
      scale: CONFIG.SCALE,
      maxPixels: 1e9
    });
    return ee.Feature(null, {
      'date': image.get('system:time_start'),
      'ET_mean': mean.get('ET'),
      'month': image.get('month'),
      'year': image.get('year')
    });
  });
  
  Export.table.toDrive({
    collection: timeSeriesData,
    description: 'Sudan_ET_Monthly_TimeSeries',
    fileFormat: 'CSV',
    folder: 'Sudan_ET_Analysis'
  });
  
  // 3. Export Trend Analysis Results
  Export.image.toDrive({
    image: trend,
    description: 'Sudan_ET_Trend_Analysis',
    scale: CONFIG.SCALE,
    region: sudan,
    maxPixels: 1e9,
    folder: 'Sudan_ET_Analysis',
    fileFormat: 'GeoTIFF'
  });
  
  // 4. Export SETI and Drought Classification
  var lastSETI = SETI_collection.select(['SETI', 'drought_class']).mosaic();
  Export.image.toDrive({
    image: lastSETI,
    description: 'Sudan_SETI_and_DroughtClass',
    scale: CONFIG.SCALE,
    region: sudan,
    maxPixels: 1e9,
    folder: 'Sudan_ET_Analysis',
    fileFormat: 'GeoTIFF'
  });
  
  // 5. Export Annual Statistics
  var annualStats = ee.List.sequence(CONFIG.START_YEAR, CONFIG.END_YEAR).map(function(year) {
    var yearCollection = ET.filter(ee.Filter.calendarRange(year, year, 'year'));
    var yearMean = yearCollection.mean();
    var stats = yearMean.reduceRegion({
      reducer: ee.Reducer.mean().combine({
        reducer2: ee.Reducer.stdDev(),
        sharedInputs: true
      }).combine({
        reducer2: ee.Reducer.minMax(),
        sharedInputs: true
      }),
      geometry: sudan,
      scale: CONFIG.SCALE,
      maxPixels: 1e9
    });
    return ee.Feature(null, {
      'year': year,
      'mean_ET': stats.get('ET_mean'),
      'stdDev_ET': stats.get('ET_stdDev'),
      'min_ET': stats.get('ET_min'),
      'max_ET': stats.get('ET_max')
    });
  });
  
  Export.table.toDrive({
    collection: ee.FeatureCollection(annualStats),
    description: 'Sudan_ET_Annual_Statistics',
    fileFormat: 'CSV',
    folder: 'Sudan_ET_Analysis'
  });
  
  // 6. Export Mean ET Raster
  Export.image.toDrive({
    image: maskedET,
    description: 'Sudan_Mean_ET_Raster',
    scale: CONFIG.SCALE,
    region: sudan,
    maxPixels: 1e9,
    folder: 'Sudan_ET_Analysis',
    fileFormat: 'GeoTIFF'
  });
};

// Map visualization
Map.centerObject(sudan, 5);

// Add layers
var meanET = ET.mean();
var maskedET = meanET.updateMask(meanET.gt(0));

Map.addLayer(
  maskedET.clip(sudan),
  etVis,
  'Mean Evapotranspiration'
);

Map.addLayer(
  ee.Image().paint(sudan, 0, 2),
  {palette: '#000000'},
  'Sudan Boundary'
);

// Add UI elements
Map.add(createLegend());
print(createTimeSeries());

// Export data
exportData();