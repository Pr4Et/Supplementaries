%Shahar added 18/5/2020
function imodModel = correct (imodModel)
imodObject = imodModel.Objects{1};
imodModel.contour =getNContours(imodObject);
dim = getMax(imodModel);
imodModel.point=dim(3);
