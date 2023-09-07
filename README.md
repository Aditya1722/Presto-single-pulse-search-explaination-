# Presto-single-pulse-search-explaination-
In the first line we assigning name to the block that if it is executed as a main prgram or module etc .

imported `hotshot` - used to profile the code block i.e statistical infos like running time and etc 

So in this block of the code we are just calling the main function as conditions under `if(0)` won`t be executed

```
def main():
    parser = OptionParser(usage)
    parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                      default=False, help="Don't make a postscript plot, just use an X-window")
    parser.add_option("-p", "--noplot", action="store_false", dest="makeplot",
                      default=True, help="Look for pulses but do not generate a plot")
    parser.add_option("-m", "--maxwidth", type="float", dest="maxwidth", default=0.0,
                      help="Set the max downsampling in sec (see below for default)")
    parser.add_option("-t", "--threshold", type="float", dest="threshold", default=5.0,
                      help="Set a different threshold SNR (default=5.0)")
    parser.add_option("-s", "--start", type="float", dest="T_start", default=0.0,
                      help="Only plot events occuring after this time (s)")
    parser.add_option("-e", "--end", type="float", dest="T_end", default=1e9,
                      help="Only plot events occuring before this time (s)")
    parser.add_option("-g", "--glob", type="string", dest="globexp", default=None,
                      help="Process the files from this glob expression")
    parser.add_option("-f", "--fast", action="store_true", dest="fast",
                      default=False, help="Use a faster method of de-trending (2x speedup)")
    parser.add_option("-b", "--nobadblocks", action="store_false", dest="badblocks",
                      default=True, help="Don't check for bad-blocks (may save strong pulses)")
    parser.add_option("-d", "--detrendlen", type="int", dest="detrendfact", default=1,
                      help="Chunksize for detrending (pow-of-2 in 1000s)")
    parser.add_option("-z", "--gzip", action="store_true", dest="gzip",
                      default=False, help="gzip the output .singlepulse files")
    (opts, args) = parser.parse_args()
```
This is the part of main fucntion in which they are defining list of options for users to use by using `optparse` python module  

End of completion we used `(opts, args) = parser.parse_args()` to tell `optparse` that we`re done

`parse_args()` returns two values:
* options, an object containing values for all of your options—e.g. if --file takes a single string argument, then options.file will be the filename supplied by the user, or None if the user did not supply that option
* args, the list of positional arguments leftover after parsing options

Here 
`Action = store`tells to take the next argument and store it to your choosen destination
and the flags true/false  just to put the dest value to true or false`
`help` - use to give some info about the command 
`dest` - destination 

```
 if len(args)==0:
        if opts.globexp==None:
            print(full_usage)
            sys.exit(0)
        else:
            args = []
            for globexp in opts.globexp.split():
                args += glob.glob(globexp)
 useffts = True
 dosearch = True
 if opts.xwin:
    pgplot_device = "/XWIN"
 else:
    pgplot_device = ""

```
In this block of code `opts.globexp` is checking the variable if given any in command line . eg ` --globexp .txt `

usefft & dosearch are assigned to be true 

if `opts.xwin` is true then pgplot_device is assignned to "/XWIN" otherwise "" 

```

    fftlen = 8192     # Should be a power-of-two for best speed
    chunklen = 8000   # Must be at least max_downfact less than fftlen
    assert(opts.detrendfact in [1,2,4,8,16,32])
    detrendlen = opts.detrendfact*1000
    if (detrendlen > chunklen):
        chunklen = detrendlen
        fftlen = int(next2_to_n(chunklen))
    blocks_per_chunk = chunklen // detrendlen
    overlap = (fftlen - chunklen) // 2
    worklen = chunklen + 2*overlap  # currently it is fftlen...

    max_downfact = 30
    default_downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300]

    if args[0].endswith(".singlepulse") or args[0].endswith(".singlepulse.gz"):
        filenmbase = args[0][:args[0].rfind(".singlepulse")]
        dosearch = False
    elif args[0].endswith(".dat"):
        filenmbase = args[0][:args[0].rfind(".dat")]
    else:
        filenmbase = args[0]
```
fftlen = 8192 (for fourier transform )

chunklen ?

* assert statement is used to check if a condition is True, and if it's not, it raises an error.
* detrendfact ? (cheking if detrendfact have the values mentioned in code )
* dtrendlen is decided
* Then dtrendlen ,chunklen,fflten all are recalculated considering dtrendfact user input (next2_to_n is a lib imported from presto )
* blocks per chunk is devision of  chunkslen and dtrendlen
* overlap , worklen is defined and max_downfact and default are constant
* checking arg
* 

```
# Don't do a search, just read results and plot
    if not dosearch:
        info, DMs, candlist, num_v_DMstr = \
              read_singlepulse_files(args, opts.threshold, opts.T_start, opts.T_end)
        orig_N, orig_dt = int(info.N), info.dt
        obstime = orig_N * orig_dt
    else:
        DMs = []
        candlist = []
        num_v_DMstr = {}

```
* In this block code we`re deciding to do a search or not by checking  if dosearch is false then it will get values of  info, DMs, candlist, num_v_DMstr from fucntion `read_singlepulse_files()` by reading single pulse files 

As they are referring to single pulse files funtion first let look into that 

```
def read_singlepulse_files(infiles, threshold, T_start, T_end):
    DMs = []
    candlist = []
    num_v_DMstr = {}
    for ii, infile in enumerate(infiles):
        if infile.endswith(".singlepulse") or infile.endswith(".singlepulse.gz"):
            filenmbase = infile[:infile.rfind(".singlepulse")]
        else:
            filenmbase = infile
        info = infodata.infodata(filenmbase+".inf")
        DMstr = "%.2f"%info.DM
        DMs.append(info.DM)
        num_v_DMstr[DMstr] = 0
        if ii==0:
            info0 = info
        if os.stat(infile)[6]:
            try:
                cands = np.loadtxt(infile)
                if len(cands.shape)==1:
                    cands = np.asarray([cands])
                for cand in cands:
                    if cand[2] < T_start: continue
                    if cand[2] > T_end: break
                    if cand[1] >= threshold:
                        candlist.append(candidate(*cand))
                        num_v_DMstr[DMstr] += 1
            except:  # No candidates in the file
                IndexError
    DMs.sort()
    return info0, DMs, candlist, num_v_DMstr
```

In this code block we are trying to read single pulse file(main objective ) 
245 and 246 , 249 , 255 ,

load all the info in cands (here doubt diff between cand and cand )

some condition and if that satisfies will store it as a candidate
`candiate` class is used to store that value in a candidate list ?









```
 # Loop over the input files
        for filenm in args:
            if filenm.endswith(".dat"):
                filenmbase = filenm[:filenm.rfind(".dat")]
            else:
                filenmbase = filenm
            info = infodata.infodata(filenmbase+".inf")
            DMstr = "%.2f"%info.DM
            DMs.append(info.DM)
            N, dt = int(info.N), info.dt
            obstime = N * dt
            # Choose the maximum width to search based on time instead
            # of bins.  This helps prevent increased S/N when the downsampling
            # changes as the DM gets larger.
            if opts.maxwidth > 0.0:
                downfacts = [x for x in default_downfacts if x*dt <= opts.maxwidth]
            else:
                downfacts = [x for x in default_downfacts if x <= max_downfact]
            if len(downfacts) == 0:
                downfacts = [default_downfacts[0]]
            if (filenm == args[0]):
                orig_N = N
                orig_dt = dt
                if useffts:
                    fftd_kerns = make_fftd_kerns(default_downfacts, fftlen)
            if info.breaks:
                offregions = list(zip([x[1] for x in info.onoff[:-1]],
                                 [x[0] for x in info.onoff[1:]]))

                # If last break spans to end of file, don't read it in (its just padding)
                if offregions[-1][1] == N - 1:
                    N = offregions[-1][0] + 1
``` 
* Here we search for `filname` if it ends with `.dat` file or not if yes then we use `rfind` to find the `.dat` file index and try to save only name of the file in `filenmbase` if not then `filename` is saved as it is in `filenmbase`

* Filenmbase we got from previous step is added with an extension `.inf` to extract the informatiom using `info.info()` and save it in `info`
* Then we extract dm info from info.dm and first add it in DMstr and then append DM values in `DMs`
* Then were extracting the info of `dt`(time step size) and `N`(no. of time bing )so we can get observation time(`obstime`)
* `maxwidth` is either given by you or default
* it checks if `maxwidth`(which you entered ) is > than 0 
  > if yes then compare it with each (`down_sample list` values)*`dt`  and if it is < than `maxwidth` than add it to the list of `downnfacts`
  > if no then store all the values which are less than the `max_downfact` which is already defined as 30
* if `downsample` list is empty then the first value of `default_sample` if filled
* This part runs only for first file it will take informtion of `N` and `dt` and store it in some variable  `orig_N` & `orig_dt`
* if `useffts` is true then it runs a function `make_fftd_kerns` by sending argument - `default_downfacts,fftlen`
As we encontered this function `make_fft_kerns` so we will discuss that before moving forward to next code part d
```
def make_fftd_kerns(downfacts, fftlen):
    fftd_kerns = []
    for downfact in downfacts:
        kern = np.zeros(fftlen, dtype=np.float32)
        # These offsets produce kernels that give results
        # equal to scipy.signal.convolve
        if downfact % 2:  # Odd number
            kern[:downfact//2+1] += 1.0
            kern[-(downfact//2):] += 1.0
        else:             # Even number
            kern[:downfact//2+1] += 1.0
            if (downfact > 2):
                kern[-(downfact//2-1):] += 1.0
        # The following normalization preserves the
        # RMS=1 characteristic of the data
        fftd_kerns.append(rfft(kern / np.sqrt(downfact), -1))
    return fftd_kerns
```
* Int his block assigning `fftd_kerns` an empty list 
* `For` loop over downfact values and creating an array of 0 of length `fftlen`(value given in main fucntion ) for each `downfact` value .
* making a shape of kernel for each downfact values depedning on if it is odd or even
* but why is it depending on 2 for even  ?
* how is it making a shape when whole fftlen is made 1
* after this we are normalising the thing
* it will always return the last kernel value no ?
*   
