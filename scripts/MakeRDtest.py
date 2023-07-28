import os
import numpy as np
from PIL import Image
import PIL
from PIL import ImageFont
from PIL import ImageDraw
import argparse
import cv2
import img2pdf

#Image helper function
#stack two or more images vertically
def vconcat_resize(img_list, interpolation
                   = cv2.INTER_CUBIC):
      # take minimum width
    w_min = min(img.shape[1]
                for img in img_list)

    # resizing images
    im_list_resize = [cv2.resize(img,
                      (w_min, int(img.shape[0] * w_min / img.shape[1])),
                                 interpolation = interpolation)
                      for img in img_list]
    # return final image
    return cv2.vconcat(im_list_resize)

# combine two images side by side
def hconcat_resize(img_list,
                   interpolation
                   = cv2.INTER_CUBIC):
    # take minimum hights
    h_min = min(img.shape[0]
                for img in img_list)

    # image resizing
    im_list_resize = [cv2.resize(img,
                       (int(img.shape[1] * h_min / img.shape[0]),
                        h_min), interpolation
                                 = interpolation)
                      for img in img_list]

    # return final image
    return cv2.hconcat(im_list_resize)

def words(STR1,STR2,outfile,n=100):
    #font = ImageFont.truetype("arial.ttf", 70)
    font = ImageFont.load_default()
    img = Image.new("RGB", (1800,300), (255,255,255))
    draw = ImageDraw.Draw(img)
    draw.text((n,10), STR1, (0,0,0),font=font)
    draw = ImageDraw.Draw(img)
    draw.text((n,150), STR2, (0,0,0),font=font)
    draw = ImageDraw.Draw(img)
    img.save(outfile)

class Variant():
  def __init__(self,chr,start,end,name,type,samples,varname,prefix,family_id):
    self.chr=chr
    self.coord=str(chr)+":"+str(start)+"-"+str(end)
    self.start=start
    self.end=end
    self.name=name
    self.type=type
    self.prefix=prefix
    self.varname=varname
    self.sample=samples
    self.samples=samples.split(",")
    self.family_id=family_id
  def pesrplotname(self,dir):
    if os.path.isfile(dir+self.family_id+"_"+self.name+".png"):
      return dir+self.family_id+"_"+self.name+".png"
    elif os.path.isfile(dir+self.family_id+"_"+self.name+".left.png") and os.path.isfile(dir+self.family_id+"_"+self.name+".right.png"):
      left = cv2.imread(dir+self.family_id+"_"+self.name+".left.png")
      right = cv2.imread(dir+self.family_id+"_"+self.name+".right.png")
      horizontal_combined = hconcat_resize([left,right])
      cv2.imwrite(dir+self.family_id+"_"+self.name+".png", horizontal_combined)
      return dir+self.family_id+"_"+self.name+".png"
    else:
      return 'Error'
  def rdplotname(self,dir,maxcutoff=float("inf")):
    if int(self.end)-int(self.start)>maxcutoff:
      medium=(int(self.end)+int(self.start))/2
      newstart=str(round(medium-maxcutoff/2))
      newend=str(round(medium+maxcutoff/2))
    else:
      newstart=self.start
      newend=self.end
    if os.path.isfile(dir+self.chr+"_"+newstart+"_"+newend+"_"+self.samples[0]+"_"+self.name+"_"+self.prefix+"_"+self.samples[0]+".jpg"):
      return dir+self.chr+"_"+newstart+"_"+newend+"_"+self.samples[0]+"_"+self.name+"_"+self.prefix+"_"+self.samples[0]+".jpg"
    elif os.path.isfile(dir+self.chr+"_"+newstart+"_"+newend+"_"+self.samples[0]+"_"+self.name+"_"+self.prefix+"_"+self.family_id+".jpg"):
      return dir+self.chr+"_"+newstart+"_"+newend+"_"+self.samples[0]+"_"+self.name+"_"+self.prefix+"_"+self.family_id+".jpg"
    else:
      return 'Error'
  def makeplot(self,pedir,rddir,outdir,flank,build="hg38"):
    if ((self.type!="INS") | (self.type!="snv") | (self.type!="indel") | (self.type!="INS:ME:SVA") |  (self.type!="INS:ME:LINE1") | (self.type!="INS:ME:ALU")):
      if int(self.end)-int(self.start)<2000:
          STR2=self.varname+" "+str(int(self.end)-int(self.start))+'bp'
      else:
          STR2=self.varname+" "+str(int((int(self.end)-int(self.start))/1000))+'kb'
    else:
      STR2=self.varname
    # get name of igv plot
    pesrplot=self.pesrplotname(pedir)
    # get name of rd plot
    rdplot=self.rdplotname(rddir, flank)
    if pesrplot!='Error' and rdplot!='Error':
        # vstack them
        igv = cv2.imread(pesrplot)
        # resize the IGV image which has white space at the bottom
        y1=0
        x1=0
        h1=3000
        w1=800
        resized_igv = igv[x1:w1, y1:h1] # get rid of white space at bottom of igv plot
        img = Image.open(rdplot) # rd plot
        img2 = img.crop((0, 230, img.size[0], img.size[1])) # crop out original RD plot annotations
        img2.save("croprd.jpg") # Harolds cropping
        rd = cv2.imread("croprd.jpg") # read it in cv2 for stacking command
        # get new annotation
        STR1=self.chr+":"+'{0:,}'.format(int(self.start))+'-'+'{0:,}'.format(int(self.end))+" (+"+build+")"
        outfile='info.jpg'
        words(STR1,STR2,outfile,100) # new Rd plot
        #img_v_resize = vconcat_resize([resized_igv,rd]) # combine rd pe and sr together
        #img_v_resize = vconcat_resize([igv,rd])
        cv2.imwrite(outdir+"rd.png", rd)
        cv2.imwrite(outdir+"igv.png", igv)
        with open(outdir+self.varname+"_denovo.png","wb") as f:
            f.write(img2pdf.convert("igv.png", "rd.png"))
    elif pesrplot!='Error' and rdplot=='Error':
        igv = cv2.imread(pesrplot)
        STR1=self.chr+":"+'{0:,}'.format(int(self.start))+'-'+'{0:,}'.format(int(self.end))+" (hg38)"
        outfile='info.jpg'
        words(STR1,STR2,outfile,50)
        #cv2.imwrite(outdir+self.varname+"_denovo.png", igv)
        with open(outdir+self.varname+"_denovo.png","wb") as f:
            f.write(img2pdf.convert("igv.png"))

class VariantInfo():
  def __init__(self,pedfile,prefix):
    self.pedfile=pedfile
    self.prefixdir={}
    if os.path.isfile(prefix):
      self.prefixfile=prefix
      self.prefix=set([])
      with open(self.prefixfile,"r") as f:
        for line in f:
          if "#" not in line:
            prefix,sample=line.rstrip().split()
            self.prefixdir[sample]=prefix
            self.prefix.add(prefix)
    else:
      self.prefix=prefix
    famdct={}
    reversedct={}
    reversedctfam={}
    self.samplelist=[]
    with open(pedfile,"r") as f:
      for line in f:
        dat=line.split()
        [fam,sample,father,mother]=dat[0:4]
        if father+","+mother not in famdct.keys():
          famdct[father+","+mother]=[sample]
        else:
          famdct[father+","+mother].append(sample)
        reversedct[sample]=father+","+mother
        reversedctfam[sample]=fam
        self.samplelist.append(sample)
    self.famdct=famdct
    self.reversedct=reversedct
    self.reversedctfam=reversedctfam

  def getprefix(self,sample):
    if self.prefixdir=={}:
      return self.prefix
    else:
      return self.prefixdir[sample]
  def getnuclear(self,sample):
    parents=self.reversedct[sample]
    if parents!="0,0":
      kids=self.famdct[parents].copy()
      kids.remove(sample)
      return sample+','+parents
    else:
      return sample

class GetVariants():
  def __init__(self,inputfile,pedfile,prefix):
    self.inputfile=inputfile
    self.variants=[]
    self.variantinfo=VariantInfo(pedfile,prefix)
    with open(inputfile,"r") as f:
      for line in f:
        if "#" not in line:
          dat=line.rstrip().split("\t")
          [chr,start,end,name,type,samples]=dat[0:6]
          sample=samples.split(',')[0]
          varname=samples.split(',')[0]+'_'+name
          if "," in sample:
            raise Exception("should only have 1 sample per variant")
          prefix=self.variantinfo.getprefix(sample)
          nuclearfam=self.variantinfo.getnuclear(sample)
          variant=Variant(chr,start,end,name,type,nuclearfam,varname,prefix,self.variantinfo.reversedctfam[sample])
          self.variants.append(variant)
  def GetRdfiles(self):
    with open(self.inputfile+".igv","w") as g:
      if self.variantinfo.prefixdir!={}:
        for prefix in self.variantinfo.prefix:
          open(self.inputfile+'_'+prefix+".txt", 'w').close()
      else:
        open(self.inputfile+'_'+self.variantinfo.prefix+".txt", 'w').close()
      for variant in self.variants:
        f=open(self.inputfile+'_'+variant.prefix+".txt",'a')
        f.write("\t".join([variant.chr,variant.start,variant.end,variant.name,variant.type,variant.sample])+'\n')
        g.write("\t".join([variant.chr,variant.start,variant.end,variant.name,variant.type,variant.sample,variant.varname])+'\n')
        f.close()

class GetDenovoPlots():
  def __init__(self,inputfile,pedfile,prefix,pedir,rddir,outdir,flank,build="hg38",GetVariantFunc=GetVariants):
    self.variants=GetVariantFunc(inputfile,pedfile,prefix).variants
    if pedir[-1]=="/":
      self.pedir=pedir
    else:
       self.pedir=pedir+"/"
    if rddir[-1]=="/":
      self.rddir=rddir
    else:
       self.rddir=rddir+"/"
    if outdir[-1]=="/":
      self.outdir=outdir
    else:
       self.outdir=outdir+"/"
    self.build=build
    self.flank=flank
  def getplots(self):
    for variant in self.variants:
      variant.makeplot(self.pedir,self.rddir,self.outdir, self.flank ,self.build)

#Main block
def main():
  parser = argparse.ArgumentParser(
      description=__doc__,
      formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('varfile')
  parser.add_argument('pedfile')
  parser.add_argument('prefix')
  parser.add_argument('flank')
  parser.add_argument('pedir')
  parser.add_argument('rddir')
  parser.add_argument('outdir')
  args = parser.parse_args()
  obj=GetDenovoPlots(args.varfile,args.pedfile,args.prefix,args.pedir,args.rddir,args.outdir,int(args.flank),"hg38",GetVariants)
  obj.getplots()
if __name__ == '__main__':
    main()