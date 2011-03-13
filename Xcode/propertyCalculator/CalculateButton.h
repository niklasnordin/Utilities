//
//  CalculateButton.h
//
//  Created by Niklas Nordin on 2011-01-03.
//  Copyright 2011 nequam. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface CalculateButton : NSObject {
    IBOutlet NSTextField *input;
    IBOutlet NSTextView *output;
}
- (IBAction)now:(id)sender;
@end
