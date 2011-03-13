//
//  propertyCalculatorAppDelegate.h
//  propertyCalculator
//
//  Created by Niklas Nordin on 2010-12-30.
//  Copyright 2010 nequam. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface propertyCalculatorAppDelegate : NSObject <NSApplicationDelegate> {
    NSWindow *window;
}

@property (assign) IBOutlet NSWindow *window;

@end
