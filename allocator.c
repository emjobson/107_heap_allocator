/*
 * File: allocator.c
 * Author: Rodrigo Grabowsky and Elliot Jobson
 * ----------------------
 * Uses segregated fit implementation, with deferred coalescing (except in myrealloc,
 * which uses immediate coalescing) and selective splitting. More details in
 * readme.txt
 *
 * Points about HONOR CODE:
 * Documentation notes from man pages used to help with writing comments.
 * Some functions adapted from macros in Section 9.9 of B&O.
 */

#include <stdlib.h>
#include <string.h>
#include "allocator.h"
#include "segment.h"
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// Heap blocks are required to be aligned to 8-byte boundary
#define ALIGNMENT 8
#define WSIZE 4 // Word and header/footer size
#define WSIZE_IN_BITS 32
#define DSIZE 8 // Double word size
#define MINQLSZ 24 // Minimum quick list block size
#define NQLS 14 // Number of quick lists
#define MAXQLSZ (MINQLSZ + (NQLS - 1)*ALIGNMENT)  // Maximum quick list block size (it is 128 for now)
#define NMISCLS 6 // Number of miscellaneous lists
#define NLISTS (NQLS + NMISCLS) // Total number of lists (or buckets): equals NQLS + NMISCLS
#define MAXSZ (1 << 30) // Maximum block size (1GB)
#define USAGE_RATIO .53 // Will split if the ratio of desired_blocksz to given_blocksz is less than this

// Global variables
static void **free_arr;
static void *heap_start;
static unsigned int mem_remaining; // Tracks amount of memory left in current heap segment (in terms of block sizes) 
static unsigned int nalloced;
static unsigned int nfreed;
static unsigned int MIN_MISC_LIST_SZ_EXP = 0; // Change this implementation later to more efficient, if time permitting
// didn't let me use const keyword

// Inline function to pack a size and alloc bits into a word
static inline unsigned int pack(unsigned int size, unsigned int alloc) {
    return (size | alloc);
}

// reads and writes  a word at address hfp*
// WARNING: careful NOT to confuse the input to these functions and pass
// a pointer to a payload instead
static inline unsigned int get(void *hfp) 
{
    return *(unsigned int *)hfp; // HFP stands for header or footer pointer
}

static inline void put(void *hfp, unsigned int val)
{
    *(unsigned int *)hfp = val;
}

// reads size field from address p*
static inline unsigned int get_size(void* hfp)
{
    return get(hfp) & ~0x7; //uses 0x7 because the size must be a multiple of 8 (so the last 3 bits won't be used)
}

// reads alloc bit from address p* ... address p must be either the header or footer
static inline bool is_alloced(void* hfp)
{
    return (get(hfp) & 0x1) == 1; //returns true if alloc bit is 1 (block is alloc'ed)
}

/* Given block ptr, compute address of its header and footer */
static inline void *hdrp(void *blockptr)
{
    return (char *)blockptr - WSIZE;
}

static inline void *ftrp(void *blockptr)
{
    return (char *)blockptr + get_size(hdrp(blockptr)) - DSIZE;
}

/* Given a pointer to the header of a block, we want to return
 * a pointer to the payload of the block (if it is alloc'd) or 
 * a pointer to the 'field' "next pointer" if it's a free block
 */
static inline void *hdr_to_payload(void *hdrp)
{
    return (char *)hdrp + WSIZE;
}

static inline void *next_to_prevptr(void *next)
{
    return (char *)next + sizeof(void *);
}


// Very efficient bitwise round of sz up to nearest multiple of mult
// does this by adding mult-1 to sz, then masking off the
// the bottom bits to compute least multiple of mult that is
// greater/equal than sz, this value is returned
// NOTE: mult has to be power of 2 for the bitwise trick to work!
static inline size_t roundup(size_t sz, int mult)
{
    return (sz + mult-1) & ~(mult-1);
}


static void setptr(void *p, void *val)
{
    *(void **)p = val;
}

// This new function dereferences the current blockp to obtain the address of the next free block
static inline void* next_free_block(void* blockp) 
{
    return *(void**)blockp;
}

// This new function takes in input of void* blockp (or void* next) to obtain the address of the prev free block
static inline void* prev_free_block(void* blockp)
{
    void* prevp = next_to_prevptr(blockp);
    return *(void**)prevp;
}

/* Given block ptr bp, compute address of next and previous blocks (adjacent in memory) */
static inline void *next_blkp(void *bp)
{
    return (char *)bp + get_size((char *)bp - WSIZE);
}

static inline void *prev_blkp(void *bp)
{
    return (char *)bp - get_size((char *)bp - DSIZE);
}



// Initialize an array of void*s of size NLISTS, with each element set as NULL
static void *initialize_array()
{
    for (int i = 0; i < NLISTS; i++) {
        free_arr[i] = NULL;
    }
    return (char*)free_arr + NLISTS * sizeof(void*); // returns address right after the end of the array
}

/* Epilogue header */
void make_epilogue()
{    
    put((char *)heap_start +  heap_segment_size() - WSIZE, pack(0, 1)); 
}

// Creates a block at the end of the allocated heap segment to mark the end of the segment
// and  a padding block (1 word) plus a prologue block (2 words) that marks the beginning of 
// the heap segment. This allows us to not account for special cases (e.g.:  running off the 
// end of the segment) when using the coalesce_block function. This idea was incorporated
// from the simple textbook implementation of the heap allocator.
void *init_prologue_epilogue(void *mem_brk)
{
    put(mem_brk, 0); /* Alignment padding */
    put((char *)mem_brk + WSIZE, pack(DSIZE, 1)); /* Prologue header */
    put((char *)mem_brk + 2*WSIZE, pack(DSIZE, 1)); /* Prologue footer */
    make_epilogue();
    return (char *)mem_brk + (3*WSIZE);
}


// This function takes in the size of a block in bytes, assuming that a valid size is entered
// (greater than MINQLSZ and under MAXSZ), and returns the index of the free list (in the free
// array) in which that block should be located, if it were a free block.
static int blocksz_to_index(unsigned int blocksz)
{
    if (blocksz <= MAXQLSZ) return (blocksz - MINQLSZ) / ALIGNMENT; // uses ALIGNMENT b/c the size for each quick list bucket increases by 8
    blocksz--;
    blocksz >>= MIN_MISC_LIST_SZ_EXP;
    int index = NQLS;
    while ((blocksz >>= 1) && (index < NLISTS - 2)) {
        index++;
    }
    return index;
}

// inserts a new free block as the first linked list element in the free_array
static void insert_freeblock(void *blockp)
{
    int size = get_size(hdrp(blockp));
    int index = blocksz_to_index(size);
    void* first = free_arr[index];
    if (first) { // if there is a free block in the list, set its prev ptr to blockp
        setptr(next_to_prevptr(first), blockp);
    }
    setptr(blockp, first); //sets the next pointer

    setptr(next_to_prevptr(blockp), NULL); //sets the prev pointer to NULL (doesn't point back to the array)
    free_arr[index] = blockp; //sets the array at that index to store a pointer to the free block

}


static void leftover_to_free_arr(void* blockp, unsigned int blocksz)
{
    put(hdrp(blockp), pack(blocksz, 0)); // set size and alloc/free status of first free block's header
    put(ftrp(blockp), pack(blocksz, 0)); //same for footer
    insert_freeblock(blockp);
}

/* The responsibility of the myinit function is to configure a new
 * empty heap. Typically this function will initialize the
 * segment (you decide the initial number pages to set aside, can be
 * zero if you intend to defer until first request) and set up the
 * global variables for the empty, ready-to-go state. The myinit
 * function is called once at program start, before any allocation 
 * requests are made. It may also be called later to wipe out the current
 * heap contents and start over fresh. This "reset" option is specifically
 * needed by the test harness to run a sequence of scripts, one after another,
 * without restarting program from scratch.
 */
bool myinit()
{
    // Set up global constant: exponent of the sizes in the range of our first misc list (powers of 2)
    // Ex.: if first misc list range is 129-256, this exp should be 7
    if (MIN_MISC_LIST_SZ_EXP == 0) MIN_MISC_LIST_SZ_EXP = (int)(log(MAXQLSZ)/log(2)); // Change this implementation later to more efficient, if time permitting

    free_arr = init_heap_segment(1); // reset heap segment to empty, no pages allocated
    heap_start = free_arr;
    if (!free_arr) return false; 
    void *mem_brk = initialize_array(); // mem_brk (ending of the currently allocated memory)  set to memory location after array
    mem_brk = init_prologue_epilogue(mem_brk); // mem_brk set to mem location right after prologue footer
    void *last_addr = (char *)heap_start + heap_segment_size();
    int size = (char *)last_addr - (char *)mem_brk - WSIZE; //calculate size of memory block to be placed in the free list
     
    void *blockp = hdr_to_payload(mem_brk); // Location of payload in alloc'd is equivalent to that of next in free'd
    
    leftover_to_free_arr(blockp, size); //initializes new free block and places it in the free array, with the leftover memory

    // Subtract the size of alignment padding, prologues header and footer, epilogue header, and freelist array
    mem_remaining = heap_segment_size() - (4 * WSIZE + NLISTS * sizeof(void *));
    return true;
}


static size_t reqsz_to_blocksz(size_t requestedsz)
{
    size_t blocksz = roundup(requestedsz + DSIZE, ALIGNMENT); //adds 8 for overhead, then rounds up to nearest multiple of 8
    if (blocksz < MINQLSZ) blocksz = MINQLSZ;
    return blocksz;   
}


/*
 *This function takes in a block size as a parameter, then searches for an appropriately-sized free
 *block by starting at the appropriate index in the free list. Iterates through each linked list,
 *continuing along the array each time a search of the current list has failed. Returns 
 *either a pointer to the payload/next field of an appropriate block, NULL if no block was found.
 */
static void* search_free_list(size_t blocksz, int* free_index)
{
    int start_index = blocksz_to_index(blocksz);
    for (int i = start_index; i < NLISTS; i++) {
        //traverse linked list
        void* curr_blockp = free_arr[i];
        while (curr_blockp) { //TODO: make sure case is handled where current bucket is already NULL
            if (get_size(hdrp(curr_blockp)) >= blocksz) {
                *free_index = i;
                return curr_blockp;
            }
            curr_blockp = next_free_block(curr_blockp);
        }
    }
    *free_index = -1; // making debugging easier
    return NULL;
}

// function that removes a free block by setting the next ptr of the prev block
// and the prev pointer of the next block
static void remove_free_block(void* blockp, int free_index)
{
    void* next_block = next_free_block(blockp);
    void* prev_block = prev_free_block(blockp); //this function takes in input of blockp, or void* next (for uniformity)

    if (!prev_block) { //if prev_block is NULL, we are removing the first element
        free_arr[free_index] = next_block;
    } else {
        setptr(prev_block, next_block); //set the next pointer of the previous block to point wherever our current next pointed
    }
    //set the previous pointer of the next block as our current block's prev pointer
    if (next_block) setptr(next_to_prevptr(next_block), prev_block);
}


// Returns a pointer to the coalesced block, and modifies the int representing the size of the coalesced block
// Returns NULL if coalescing was not possible
static void* coalesce_block(void* bp, int* coalescedsz)
{
    // Modifying the function so that it is also able to return a coalesced block that is allocated
    // instead of just blocks that are free (by setting alloc/free bit to af_bit instead of 0)
    // For use in myrealloc
    unsigned int af_bit = is_alloced(hdrp(bp)) ? 1 : 0; 
    int leftsz = 0;
    // coalesce left
    void *curr = prev_blkp(bp);
    void* new_bp = bp;
    while (!is_alloced(ftrp(curr))) { //if previous block is free
        int currsz = get_size(ftrp(curr));
        leftsz += currsz;
        remove_free_block(curr, blocksz_to_index(currsz));
        new_bp = curr;
        curr = prev_blkp(curr);
    }
    
    int rightsz = 0;
    // coalesce right
    curr = next_blkp(bp);
    while (!is_alloced(hdrp(curr))) {
        int currsz = get_size(hdrp(curr));
        rightsz += currsz;
        remove_free_block(curr, blocksz_to_index(currsz));
        curr = next_blkp(curr);
    }

    if (leftsz != 0 || rightsz != 0) { // If any coalescing has occurred
        int basesz = get_size(hdrp(bp));
        if (af_bit == 0) remove_free_block(bp, blocksz_to_index(basesz)); // Remove the original block
	int totalsz = basesz + leftsz + rightsz;
        put(hdrp(new_bp), pack(totalsz, af_bit)); // Set header size and alloc/free bit of newly allocated block
	put(ftrp(new_bp), pack(totalsz, af_bit)); // Set footer size and alloc/free bit

        *coalescedsz = totalsz;
        return new_bp;
    }
    return NULL;
}

// function that coalesces wherever possible. returns a pointer to a coalesced block of the appropriate size, if possible
// handles iteration through free_arr
static void* def_coalesce(size_t blocksz)
{
    void* free_block = NULL;
    for (int i = 0; i < NLISTS; i++) {
        void* curr = free_arr[i];
        while (curr) {
            int coalescedsz = 0;
            void* new_coalesce = coalesce_block(curr, &coalescedsz);
            if (new_coalesce) {
                if (!free_block && coalescedsz >= blocksz) {
                    free_block = new_coalesce;
                } else { //else need to place coalesced block back into free_arr
                    insert_freeblock(new_coalesce);
                }   
            }
            curr = next_free_block(curr);
         }
    }
    return free_block;
}


static void *request_memory(size_t blocksz)
{
    void *old_epilogue = (char *)heap_start + heap_segment_size() - WSIZE;
    int npages = roundup(2*blocksz, PAGE_SIZE) / PAGE_SIZE;
    // testing with EJs idea of allocating new mem as if they were going to be one more request of the same size
    extend_heap_segment(npages);
    // now, make new epilogue
    make_epilogue();
    // make allocated block
    void *ret = hdr_to_payload(old_epilogue);
    put(old_epilogue, pack(blocksz, 1)); // header of new alloc block
    put(ftrp(ret), pack(blocksz, 1)); //footer of new alloc block
    
    // place free memory as new free block in free list
    void* free_blockp = (char*)ret + blocksz;
    unsigned int leftover_free_size = npages * PAGE_SIZE - blocksz;
    mem_remaining += leftover_free_size;
    leftover_to_free_arr(free_blockp, leftover_free_size);
    
    return ret;
}

static void init_malloc_blk(void *frblkp)
{
    put(hdrp(frblkp), pack(get_size(hdrp(frblkp)), 1)); 
    put(ftrp(frblkp), pack(get_size(ftrp(frblkp)), 1));
}

// Usage ratio : 53%
static void* check_splitting(void* blockp, unsigned int wanted_blocksz, unsigned int given_blocksz)
{
    if (given_blocksz - wanted_blocksz > 24 && ((double)wanted_blocksz / (double)given_blocksz) <= USAGE_RATIO) { //if a new block can be formed, and it's worth it, split
        //update return blockp's header and footer
        put(hdrp(blockp), pack(wanted_blocksz, 1));
        put(ftrp(blockp), pack(wanted_blocksz, 1));

        //initialize leftover free block
        void* leftover_freeblk = (char*)blockp + wanted_blocksz;
        size_t leftoversz = given_blocksz - wanted_blocksz;
        put(hdrp(leftover_freeblk), pack(leftoversz, 0));
        put(ftrp(leftover_freeblk), pack(leftoversz, 0));
        
        insert_freeblock(leftover_freeblk); //place leftover_freeblk into free_arr
        mem_remaining -= wanted_blocksz; //decrease the memory remaining
        return blockp;
    } 
    //splitting not necessary, retutn NULL as signal that splitting did not occur
    return NULL;
}

static void *alloc_routine(void *bp, unsigned int blocksz) 
{
    void *split_block = check_splitting(bp, blocksz, get_size(hdrp(bp)));
    if (split_block) return split_block;
    init_malloc_blk(bp);
    mem_remaining -= get_size(hdrp(bp));
    return bp;
}

void *mymalloc(size_t requestedsz)
{
    // if size is 0 or greater than upper limit (1GB), then malloc() returns either
    // NULL, or a unique pointer value that can later be successfully passed to free().
    if (requestedsz == 0 || requestedsz > MAXSZ) return NULL;
    size_t blocksz = reqsz_to_blocksz(requestedsz);
    if (blocksz > mem_remaining) { // Need to ask for more pages
        nalloced++;
        return request_memory(blocksz);

    } else { //enough total unallocated memory exists to satisfy the request        
        //search free list for an appropriately-sized free block
        int free_index;
        void* blockp = search_free_list(blocksz, &free_index);
        
        if (blockp) { //if a free block has been found, remove it from the free_arr data structure
            remove_free_block(blockp, free_index);
            nalloced++;
            return alloc_routine(blockp, blocksz);
	} else {  //attempt to gain a large enough free block by coalescing (save a pointer to the first block created that works
            blockp = def_coalesce(blocksz);
            nalloced++;
            if (blockp) return alloc_routine(blockp, blocksz);
        }
        nalloced++;
        return request_memory(blocksz);
    }
    return NULL;
}


// The free() function frees the memory space pointed to by ptr, which must have been returned by a previous call to mymalloc() 
// or myrealloc(). Otherwise, or if free(ptr) has already been called before, undefined behavior occurs. If ptr is NULL, 
// no operation is performed.
void myfree(void *ptr)
{
    if (!ptr) return;
    
    put(hdrp(ptr), pack(get_size(hdrp(ptr)), 0));
    put(ftrp(ptr), pack(get_size(ftrp(ptr)), 0));
    insert_freeblock(ptr);
    
    nfreed++;
}

void *myrealloc(void *oldptr, size_t newsz)
{
    if (!oldptr) return mymalloc(newsz);
    if (newsz == 0) {
        myfree(oldptr);
        return NULL;
    }
    size_t newblocksz = reqsz_to_blocksz(newsz);
    size_t oldblocksz = get_size(hdrp(oldptr));

    if (newblocksz < oldblocksz) {
        void *splitptr = check_splitting(oldptr, newblocksz, oldblocksz);
        if (splitptr) return splitptr;
        return oldptr;
    }
    if(newblocksz == oldblocksz) return oldptr;
    
    // If control passes two if statements above, it means newblocksz is greater than oldblocksz
    // So, we first try coalescing to satisfy request
    int coalescedsz = 0;
    void *to_free = oldptr;
    void *newptr = coalesce_block(oldptr, &coalescedsz);
    if (newptr) {
        if (coalescedsz >= newblocksz) {
            memmove(newptr, oldptr, oldblocksz - DSIZE);
            return newptr;
        }
        to_free = newptr;
    }

    // Experiment with actually not calling malloc and free here
    newptr = mymalloc(newsz);
    memcpy(newptr, oldptr, oldblocksz - DSIZE);
    myfree(to_free);
    return newptr;
}


// validate_heap is your debugging routine to detect/report
// on problems/inconsistency within your heap data structures
bool validate_heap()
{
    bool frees_valid = true;
    bool blk_sizes_valid = true;
    bool blocks_aligned = true;
    bool valid_pointers = true;

    void *curr;
    //check that free list blocks are actually free
    for (int i = 0; i < NLISTS; i++) {
        curr = free_arr[i];
        while (curr) {
            void* next_block = next_free_block(curr);
            void* prev_block = prev_free_block(curr);
            if (is_alloced(hdrp(curr))) frees_valid = false;
            if (next_block != NULL && (char *)next_block > (char *)heap_start + heap_segment_size() - WSIZE) valid_pointers = false;
            if (prev_block != NULL && (char *)prev_block > (char *)heap_start + heap_segment_size() - WSIZE) valid_pointers = false;
            curr = next_block;
        }
    }

    //make sure that block sizes in heap are between MINQLSZ and MAXSZ
    curr = (char *)free_arr + NLISTS * sizeof(void *) + WSIZE; // go to end of array, then move past
	// do while get_size is not 0                                  // padding, prologue and a header
    while (get_size(hdrp(curr)) != 0) { // while the epilogue has not been reached
        if (get_size(hdrp(curr)) < MINQLSZ || get_size(hdrp(curr)) > MAXSZ) blk_sizes_valid = false;
        if ((size_t)curr % ALIGNMENT != 0) blocks_aligned = false;
        curr = next_blkp(curr);
    }

    if (!frees_valid) {
        printf("Validate error: frees not valid (there are allocated blocks in free list)\n");
        return false;
    }
    if (!blk_sizes_valid) {
        printf("Validate error: block sizes not valid\n");
        return false;
    }
    if (!blocks_aligned) {
        printf("Validate error: blocks not aligned\n");
        return false;
    }
    if (!valid_pointers) {
        printf("Validate error: pointers not valid\n");
        return false;
    }
    return true; 
}








