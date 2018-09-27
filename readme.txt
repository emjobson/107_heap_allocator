NAME(s) 
Elliott Jobson: emjobson
Rodrigo Grabowsky: rmgrab

DESIGN 

DATA STRUCTURES
We used a segregated list implementation, with 20 buckets in our array. The first 14 arrays point to linked lists containing free blocks of sizes 24 for index 0, 
32 for index 1, and so on, until index 13 which points to blocks of size 128 (note: for these buckets, all blocks in a list have the same size). The lists in each
of the next 6 buckets store free blocks with sizes ranging between two powers of 2. So, for example the list at index 14 stores sizes 129 (2^7+1) through 256 (2^8),
inclusive. Index 19, stores blocks of sizes 4097 (2^12) through 1GB (2^30), which is the maximum possible size for a block. Each memory block contains at least space
for a header (1 word), a pointer to the next block in the list (2 words), a previous pointer (2 words) and a footer (1 word). Even allocated blocks, which don't need
a next and previous pointer, need at least 24 bytes of memory because if they did get free'd by a client, we'd need to be able to store them in the linked lists
and pointers would need to be created). In the heap segmented we allocated, the array of void*'s is comes first, followed by a padding block (to help satisy alignment 
requirement) and, a prologue header and footer (to help avoid special cases when coalescing the heap). In the last word of the heap segment we store an epilogue block. 

ALGORITHMS
To avoid internal fragmentation, we used selective splitting. To combat false fragmentation and external fragmentation we deferred coalescing until we received an 
allocation request we could not satisfy with any of the blocks in the segregated free list. Both these algorithms also gave us the best throughputs compared to use 
other algorithms we tried. When searching for a free block in our array of free lists, we use a first fit approach. We also used first fit when coalescing the heap: 
as soon as we coalesced specific blocks into a block large enough to satisfy the current allocation request, we saved a pointer to that block. But we also continued 
coalescing all adjacent free blocks in the heap. Whenever coalescing fails to produce a block large enough to fulfill a request we request enough memory from the 
segment module so that a potential subsequent request of the same size could be satisfied without having to request more pages. 

RATIONALE 

We used segregated free lists, because it is faster to search than an implicit or explicit free list. This is because by initially using the block size to hash into
a bucket we are searching through a much smaller free list than if all free blocks were stored in a single free list. We refer to the lists in buckets 0-13 as 'quick
lists'. The segregated free list + quick list design also allows us to use the first fit approach without the inefficiencies that first fit has in an explicit free
list. First fit in a quick list gives us a perfect match if the list is not empty and first fit in our lists in general approximates best fit because the free blocks
are organized by size ranges. We might not get the best block possible in a miscellaneous list (index>13) but at worst case we will select a block that is ~2x as 
large as the size we need (unless the list we hash into is empty, in which case, a 'good' fit would not even be possible to begin with). We used selective splitting
to speed up the program (not executing splitting every time it is possible) and to balance internal and external fragmentation (allowing some padding as opposed to
many small free blocks). We considered using a tree implementation for even better throughput but thought it would be complex to write in the limited amount of time
we had for this project. We used deferred coalescing (as opposed to immediate) because we thought it'd give us better throughput, but it turns out it ALSO helps
with utilization. 

OPTIMIZATION 
Firstly, we declared constants for everything we could: there's one for # of quick lists, # of miscellaneous lists, minimum size of block we store in a misc list, etc.
In the case of this segregated list data structure, we experimented with the constants and found that our numbers were a sweet spot between utilization and throughput.
We also implemented myrealloc without simply calling myfree and myalloc when it was not necessary (the block was simply returned if it already contained enough space
and we attempt immediate coalescing of surrounding blocks before calling myalloc. We compiled the code with the best gcc optimization flag we found, -O3. For selective
splitting we defined a usage ratio constant. Whenever we are using a block that wastes (1-usage ratio) of space (proportional to the block size) we split it. We
experimented with this constant, extensively, finding that 53% was the optimal value for both throughput and utilization, to the nearest 1%. We experimented with
many different potential improvements, including, but not limited to, immediate coalescing (to see if it gave us better utilization -- it didn't) and requesting
more than 1 page of memory in myinit. Used callgrind to remove inefficiencies -- for example, it helped us notice we had to remove expensive assert statements we
had included during debugging; it improved our throughput by many percentage points.

EVALUATION 
Weaknesses:
- Due to the fact that we use a linked list implementation and use coalescing, we needed a minimum size of 24 bytes (header, footer and 2 pointers) for every
allocated and free block; this means that for very small blocks our implementation has very low utilization rates
- Has unusually low utilization in a few specific testing scripts, namely: merge-pattern (18%), reassemble-trace (31%) and coalesce-pattern (33%). Also,
some scripts lead to very low throughputs, namely: surprise-pattern and random-pattern. I don't quite know if this is common for industry-level memory
allocators, but it seems that out implementation is very uneven accross different usages 
- Didn't try all optimizations we had in mind due to lack of time. I read in page 865 of B&O that the GNU malloc package in the C standard lib uses a type
of segregated fit implementation, so there's definitely something we are overlooking that could make both our utilization and throughput higher
- A binary search tree, or red-black tree might be faster

Strengths: 
- Exceeds benchmark for throughput (in most myth machines, we got 94%, in busy ones, it was a little slower)
- Has benchmark leve utilization (in all of our final alloctests we got 65%)
- Using next and previous pointers (doubly-linked list) allows us to remove a block from the segregated free list in constant time. With a singly-linked
list we'd have to search for the data structure to find the preceding block. 
- Faster look-up of free blocks than an implicit or explicit free list would allow
- Deferred coalescing and selective splitting were the best algorithms we could find to deal with fragmentation issues
- Headers and footers encode both size and info on whether block is allocated or freed with only 4 bytes
- Function that hashes from block size to array index is fast because it takes advantage of __builtin_clz function

REFERENCES
- Used section 9.9 of B&O textbook extensively. We incorporated many of their macros such as GET and PUT into our code as one-line functions. Also had the idea
for a segregated fit implementation from pgs. 864-865 of this textbook. We then improved their suggested implementation with algorithms such as deferred coalescing
and selective splitting.
- Read part of IBM research paper "Scalability of Dynamic Storage Allocation Algorithms" and adapted some ideas from their Multiple Free List Fit 1 algorithm,
such as using around 5 miscellaneous lists (we ended up using 6). 

