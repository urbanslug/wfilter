pub struct EditCigar {
    operations: String,
    max_operations: i32,
    begin_offset: i32,
    end_offset: i32,
    score: i32,
}

impl EditCigar {
    pub fn new(pattern_length: i32, text_length: i32) -> Self {
        // TODO: will we need this?
        // let operations: String = (char*)wflambda_mm_allocator_malloc(mm_allocator,edit_cigar->max_operations);
        let operations = String::new();
        let max_operations: i32 = pattern_length + text_length;
        let begin_offset: i32 = max_operations - 1;
        let end_offset: i32 = max_operations;
        let score: i32 = i32::MIN;

        Self {
            operations,
            max_operations,
            begin_offset,
            end_offset,
            score,
        }
    }
}

pub fn edit_cigar_allocate(pattern_length: i32, text_length: i32) -> EditCigar {
    EditCigar::new(pattern_length, text_length)
}
