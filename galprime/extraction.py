import sys
import numpy as np
from photutils.isophote import Ellipse, EllipseGeometry

from photutils.morphology import data_properties


class IsophoteFitter:
    def __init__(self, data, config={}, mask=None, centre=None, **kwargs):
        self.data = data
        self.config = config
        self.mask = mask
        self.centre = centre

        self.fail_count, self.max_fails = 0, 100

        self.linear = config.get("MODEL", {}).get("LINEAR", False)
        self.step = config.get("MODEL", {}).get("STEP", 0.1)

        if self.centre == None:
            self.centre = (self.data.shape[0] / 2, self.data.shape[1] / 2)
    

    def _fit_single(self, geometry, **kwargs):
        try:
            flux = Ellipse(self.data, geometry)
            fitting_list = flux.fit_image()
        except KeyboardInterrupt:
            sys.exit(1)
        except (RuntimeError, ValueError, OverflowError, IndexError):
            self.fail_count += 1
            if self.fail_count > self.max_fails:
                raise RuntimeError("Too many failed fit attempts")

    def fit(self):
        # First try fit using geometry guess
        try:
            morph_cat = data_properties(self.data, mask=self.mask, background=None)
            r = 2.0
            pos = (morph_cat.xcentroid, morph_cat.ycentroid)
            a = morph_cat.semimajor_sigma.value * r
            b = morph_cat.semiminor_sigma.value * r
            theta = morph_cat.orientation.value

            geometry = EllipseGeometry(x0=pos[0], y0=pos[1], sma=a, 
                                       eps=1 - b / a, pa=morph_cat.orientation.value)
            
            isolist = self._fit_single(geometry)
            

        except Exception as e:
            self.fail_count += 1

        # If that fails, try a brute force search
        try:
            for angle in range(0, 180, 45):
                for sma in range(2, 26, 5):
                    for eps in (0.3, 0.5, 0.9):
                        geometry = EllipseGeometry(x0=self.centre[0], y0=self.centre[1], 
                                                   sma=sma, eps=eps, pa=angle)
                        try:
                            isolist = self._fit_single(geometry)
                            
                            if len(isolist) > 0:
                                break
                            
                        except RuntimeError as e:
                            return []
                        except Exception as e:
                            continue

        except Exception as e:
            self.fail_count += 1
        
        
